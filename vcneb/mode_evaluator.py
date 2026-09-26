"""Calculator-independent, recoverable energy/gradient evaluations for modes.

Each point receives a unique immutable directory. Geometry is checked before
any calculator call; successful outputs can be reused after a process restart
only when the reference chart, declared calculator contract and pressure match.
This module does not schedule jobs or change a calculator's physical settings.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Callable

import numpy as np
from ase import Atoms
from ase.io import write

from .reference_cell import ReferenceCellCoordinates


Array = np.ndarray
CalculatorFactory = Callable[[int, Atoms, Path], object]
ResultValidator = Callable[[int, Atoms, Path], dict]


def _fingerprint(value: dict) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()


def _write_json(path: Path, value: dict) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    temporary.replace(path)


class CalculatorModeEvaluator:
    """Callable ``coordinates -> (H, dH/dcoordinates)`` through ASE.

    A fresh factory-created calculator is used for each unique coordinate
    vector. ``calculator_id`` must identify the backend and its immutable
    functional/pseudo/basis/k-point/SCF settings; this class records but cannot
    independently prove that user-supplied identity. Existing workdirs are
    refused unless ``resume=True`` with an exact contract match. Incomplete
    evaluation directories are preserved and retried in *new* directories.
    """

    def __init__(
        self,
        chart: ReferenceCellCoordinates,
        calculator_factory: CalculatorFactory,
        *,
        workdir: str | Path,
        calculator_id: str,
        pressure_eV_per_A3: float = 0.0,
        minimum_distance_A: float = 0.0,
        result_validator: ResultValidator | None = None,
        validation_id: str | None = None,
        resume: bool = False,
    ) -> None:
        if not isinstance(chart, ReferenceCellCoordinates) or not callable(calculator_factory):
            raise TypeError("chart and calculator_factory are required")
        if not calculator_id.strip():
            raise ValueError("calculator_id must identify the fixed physical settings")
        if (result_validator is None) != (validation_id is None):
            raise ValueError("result_validator and a nonempty validation_id must be supplied together")
        if result_validator is not None and (not callable(result_validator) or not validation_id.strip()):
            raise ValueError("result_validator must be callable and validation_id nonempty")
        pressure = float(pressure_eV_per_A3)
        minimum_distance = float(minimum_distance_A)
        if not np.isfinite(pressure) or not np.isfinite(minimum_distance) or minimum_distance < 0.0:
            raise ValueError("pressure and minimum distance must be finite; distance must be nonnegative")
        reference = chart.reference_atoms
        contract = {
            "schema_version": 1,
            "coordinate_convention": "reference_cell_atomic_displacements_A_plus_symmetric_strain_voigt",
            "reference": {
                "numbers": reference.get_atomic_numbers().tolist(),
                "cell_A": reference.cell.array.tolist(),
                "scaled_positions": reference.get_scaled_positions(wrap=False).tolist(),
                "pbc": reference.pbc.tolist(),
            },
            "calculator_id": calculator_id,
            "validation_id": validation_id,
            "pressure_eV_per_A3": pressure,
            "minimum_distance_A": minimum_distance,
        }
        contract["contract_sha256"] = _fingerprint(contract)
        self.chart = chart
        self.factory = calculator_factory
        self.result_validator = result_validator
        self.workdir = Path(workdir)
        self.contract = contract
        self.pressure = pressure
        self.minimum_distance = minimum_distance
        self._cache: dict[str, tuple[float, Array]] = {}
        self._next_index = 0
        self.n_new_evaluations = 0
        if self.workdir.exists():
            if not resume:
                raise FileExistsError(f"mode evaluation workdir already exists: {self.workdir}")
            recorded = json.loads((self.workdir / "contract.json").read_text(encoding="utf-8"))
            if recorded != contract:
                raise ValueError("mode evaluation resume contract changed")
            for directory in sorted(self.workdir.glob("eval-*-*")):
                if not directory.is_dir():
                    continue
                try:
                    index = int(directory.name.split("-", 2)[1])
                except (IndexError, ValueError) as exc:
                    raise ValueError(f"malformed mode evaluation directory: {directory}") from exc
                self._next_index = max(self._next_index, index + 1)
                result_path = directory / "result.json"
                if not result_path.is_file():
                    continue  # Preserve partial/failed output, retry in a new dir.
                result = json.loads(result_path.read_text(encoding="utf-8"))
                coordinate = np.asarray(result["coordinates"], dtype=float)
                digest = hashlib.sha256(coordinate.tobytes()).hexdigest()
                gradient = np.asarray(result["gradient"], dtype=float)
                if (result.get("contract_sha256") != contract["contract_sha256"]
                        or result.get("coordinate_sha256") != digest
                        or coordinate.shape != (chart.coordinate_count,)
                        or gradient.shape != coordinate.shape
                        or (validation_id is not None and not isinstance(result.get("validation"), dict))
                        or not np.all(np.isfinite(coordinate))
                        or not np.all(np.isfinite(gradient))
                        or not np.isfinite(result["enthalpy_eV"])):
                    raise ValueError(f"corrupt cached mode evaluation: {result_path}")
                if digest in self._cache:
                    raise ValueError(f"duplicate completed mode evaluation for coordinate hash {digest}")
                self._cache[digest] = float(result["enthalpy_eV"]), gradient
        else:
            self.workdir.mkdir(parents=True)
            _write_json(self.workdir / "contract.json", contract)

    def __call__(self, coordinates: Array) -> tuple[float, Array]:
        values = np.asarray(coordinates, dtype=float).reshape(-1)
        atoms = self.chart.to_atoms(values)  # validates size, finiteness and cell SPD
        distances = atoms.get_all_distances(mic=True)
        np.fill_diagonal(distances, np.inf)
        closest = float(np.min(distances))
        if closest < self.minimum_distance:
            raise ValueError(
                f"mode trial minimum atomic distance {closest:.6f} A is below "
                f"{self.minimum_distance:.6f} A; calculator was not called"
            )
        digest = hashlib.sha256(values.tobytes()).hexdigest()
        if digest in self._cache:
            energy, gradient = self._cache[digest]
            return energy, gradient.copy()
        index = self._next_index
        self._next_index += 1
        directory = self.workdir / f"eval-{index:06d}-{digest[:12]}"
        directory.mkdir()
        write(directory / "geometry.extxyz", atoms, format="extxyz")
        try:
            atoms.calc = self.factory(index, atoms, directory)
            energy, gradient = self.chart.energy_and_gradient(
                values, atoms, pressure_eV_per_A3=self.pressure,
            )
            forces = np.asarray(atoms.get_forces(), dtype=float)
            stress = np.asarray(atoms.get_stress(voigt=False), dtype=float)
            validation = (None if self.result_validator is None
                          else self.result_validator(index, atoms, directory))
            if validation is not None and not isinstance(validation, dict):
                raise ValueError("result_validator must return a JSON-object-like dictionary")
            result = {
                "status": "complete", "contract_sha256": self.contract["contract_sha256"],
                "coordinate_sha256": digest, "coordinates": values.tolist(),
                "enthalpy_eV": energy, "gradient": gradient.tolist(),
                "volume_A3": float(atoms.get_volume()),
                "maximum_atomic_force_eV_per_A": float(np.max(np.linalg.norm(forces, axis=1))),
                "stress_eV_per_A3": stress.tolist(),
                "minimum_distance_A": closest,
                "validation": validation,
            }
            _write_json(directory / "result.json", result)
        except Exception as exc:
            _write_json(directory / "failure.json", {
                "status": "failed", "coordinate_sha256": digest,
                "contract_sha256": self.contract["contract_sha256"],
                "exception_type": type(exc).__name__, "message": str(exc),
            })
            raise
        self._cache[digest] = energy, gradient.copy()
        self.n_new_evaluations += 1
        return energy, gradient.copy()


__all__ = ["CalculatorModeEvaluator"]
