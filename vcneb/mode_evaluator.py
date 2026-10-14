"""Calculator-independent, recoverable energy/gradient evaluations for modes.

Each point receives a unique immutable directory. Geometry is checked before
any calculator call; successful outputs can be reused after a process restart
only when the reference chart, declared calculator contract and pressure match.
This module does not schedule jobs or change a calculator's physical settings.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from typing import Callable

import numpy as np
from ase import Atoms
from ase.io import write

from .reference_cell import ReferenceCellCoordinates
from .clamped_mode_coordinates import ClampedModeCoordinates
from .biaxial_mode_coordinates import BiaxialClampedModeCoordinates


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
    Optional ``max_evaluation_attempts`` bounds cumulative reserved attempts,
    including failures, independently of the physical cache contract. It is
    not a raw-SCF count or a distributed scheduler/resource budget.
    """

    def __init__(
        self,
        chart: ReferenceCellCoordinates | ClampedModeCoordinates | BiaxialClampedModeCoordinates,
        calculator_factory: CalculatorFactory,
        *,
        workdir: str | Path,
        calculator_id: str,
        pressure_eV_per_A3: float = 0.0,
        minimum_distance_A: float = 0.0,
        result_validator: ResultValidator | None = None,
        validation_id: str | None = None,
        resume: bool = False,
        max_evaluation_attempts: int | None = None,
    ) -> None:
        if not isinstance(chart, (ReferenceCellCoordinates, ClampedModeCoordinates, BiaxialClampedModeCoordinates)) or not callable(calculator_factory):
            raise TypeError("chart and calculator_factory are required")
        if not calculator_id.strip():
            raise ValueError("calculator_id must identify the fixed physical settings")
        if max_evaluation_attempts is not None and (
            isinstance(max_evaluation_attempts, (bool, np.bool_))
            or not isinstance(max_evaluation_attempts, (int, np.integer))
            or max_evaluation_attempts < 0
        ):
            raise ValueError("max_evaluation_attempts must be a nonnegative integer or None")
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
        if isinstance(chart, (ClampedModeCoordinates, BiaxialClampedModeCoordinates)):
            # Preserve the existing symmetric-chart resume contract exactly.
            # A clamped run instead freezes its actual active cell subspace.
            contract["coordinate_convention"] = chart.coordinate_convention
            contract["mechanical_boundary"] = chart.mechanical_contract
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
        self._attempts = 0
        self._budget = (None if max_evaluation_attempts is None else {
            "schema_version": 1,
            "max_evaluation_attempts": int(max_evaluation_attempts),
            "physical_contract_sha256": contract["contract_sha256"],
        })
        self._reservation_digests: dict[int, str] = {}
        if self.workdir.exists():
            if not resume:
                raise FileExistsError(f"mode evaluation workdir already exists: {self.workdir}")
            recorded = json.loads((self.workdir / "contract.json").read_text(encoding="utf-8"))
            if recorded != contract:
                raise ValueError("mode evaluation resume contract changed")
            self._load_attempt_reservations()
            for directory in sorted(self.workdir.glob("eval-*-*")):
                if not directory.is_dir():
                    continue
                try:
                    index = int(directory.name.split("-", 2)[1])
                except (IndexError, ValueError) as exc:
                    raise ValueError(f"malformed mode evaluation directory: {directory}") from exc
                if self._budget is not None:
                    digest = self._reservation_digests.get(index)
                    if digest is None or directory.name != f"eval-{index:06d}-{digest[:12]}":
                        raise ValueError("mode evaluation directory lacks its budget reservation")
                else:
                    self._attempts += 1
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
                        or (self._budget is not None and digest != self._reservation_digests[index])
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
            self._next_index = max(self._next_index, self._attempts)
        else:
            self.workdir.mkdir(parents=True)
            _write_json(self.workdir / "contract.json", contract)
            if self._budget is not None:
                _write_json(self.workdir / "evaluation-budget.json", self._budget)

    @property
    def fixed_coordinate_covectors(self) -> np.ndarray:
        """Prescribed controls; conditional release must annihilate these columns."""
        if isinstance(self.chart, BiaxialClampedModeCoordinates):
            return self.chart.fixed_coordinate_covectors
        return np.empty((self.chart.coordinate_count, 0))

    @property
    def n_evaluation_attempts(self) -> int:
        """Cumulative reservations, including failures; not raw SCF calls.

        In capped namespaces a reservation interrupted before the factory is
        conservatively consumed. Legacy uncapped namespaces count preserved
        evaluation directories. Completed cache hits consume no reservation.
        """
        return self._attempts

    def _load_attempt_reservations(self) -> None:
        path = self.workdir / "evaluation-budget.json"
        markers = tuple(self.workdir.glob("attempt-*.json"))
        if self._budget is None:
            if path.exists() or markers:
                raise ValueError("mode evaluation attempt budget changed or removed")
            return
        if (not path.is_file()
                or _fingerprint(json.loads(path.read_text(encoding="utf-8"))) != _fingerprint(self._budget)):
            raise ValueError("mode evaluation attempt budget changed or missing")
        for marker in markers:
            data = json.loads(marker.read_text(encoding="utf-8"))
            if not isinstance(data, dict):
                raise ValueError("corrupt mode evaluation attempt reservation")
            index = data.get("evaluation_index")
            values = np.asarray(data.get("coordinates"), dtype=float)
            digest = hashlib.sha256(values.tobytes()).hexdigest()
            if (isinstance(index, bool) or not isinstance(index, int) or index < 0
                    or marker.name != f"attempt-{index:06d}.json"
                    or type(data.get("schema_version")) is not int or data["schema_version"] != 1
                    or data.get("budget_sha256") != _fingerprint(self._budget)
                    or data.get("coordinate_sha256") != digest
                    or values.shape != (self.chart.coordinate_count,)
                    or not np.all(np.isfinite(values))):
                raise ValueError("corrupt mode evaluation attempt reservation")
            self._reservation_digests[index] = digest
        count = len(self._reservation_digests)
        if (set(self._reservation_digests) != set(range(count))
                or count > self._budget["max_evaluation_attempts"]):
            raise ValueError("mode evaluation attempt reservations exceed or omit budget slots")
        self._attempts = count

    def _reserve_attempt(self, index: int, values: Array, digest: str) -> None:
        if self._budget is None:
            return
        if self._attempts >= self._budget["max_evaluation_attempts"]:
            raise RuntimeError("mode evaluation attempt budget exhausted; calculator was not called")
        record = {
            "schema_version": 1, "evaluation_index": index,
            "coordinate_sha256": digest, "coordinates": values.tolist(),
            "budget_sha256": _fingerprint(self._budget),
        }
        # Exclusive indexed slots refuse stale/concurrent writers. Never
        # remove/refund interrupted reservations or rewrite a prior slot.
        with (self.workdir / f"attempt-{index:06d}.json").open("x", encoding="utf-8") as stream:
            json.dump(record, stream, allow_nan=False)
            stream.flush()
            os.fsync(stream.fileno())
        self._reservation_digests[index] = digest
        self._attempts += 1

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
        self._reserve_attempt(index, values, digest)
        self._next_index += 1
        directory = self.workdir / f"eval-{index:06d}-{digest[:12]}"
        directory.mkdir()
        if self._budget is None:
            self._attempts += 1
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
