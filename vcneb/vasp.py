"""Small VASP helpers for VC-NEB examples."""

from __future__ import annotations

import os
import hashlib
import json
import re
import shutil
from pathlib import Path
from typing import Mapping, Optional, Sequence

import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.singlepoint import SinglePointCalculator
from ase.calculators.vasp import Vasp
from ase.io import write

from .provenance import endpoint_structure_record


REQUIRED_VCNEB_STATIC_PARAMETERS = {
    "ibrion": -1,
    "nsw": 0,
    "isif": 2,
    "isym": 0,
}


def potcar_dataset_labels(path: str | Path) -> list[str]:
    """Return dataset labels (for example ``Ba_sv``) from a POTCAR."""

    text = Path(path).read_text(encoding="latin-1")
    labels = []
    for title in re.findall(r"^\s*TITEL\s*=\s*(.+)$", text, flags=re.MULTILINE):
        candidates = re.findall(r"\b[A-Z][a-z]?(?:_[A-Za-z0-9]+)?\b", title)
        candidates = [item for item in candidates if not item.startswith("PAW")]
        if not candidates:
            raise ValueError(f"cannot determine POTCAR dataset from TITEL={title!r}")
        labels.append(candidates[0])
    if not labels:
        raise ValueError(f"no TITEL records found in POTCAR: {path}")
    return labels


def potcar_setups(path: str | Path) -> dict[str, str]:
    """Translate explicit POTCAR labels into ASE setup suffixes."""

    setups: dict[str, str] = {}
    for label in potcar_dataset_labels(path):
        match = re.fullmatch(r"([A-Z][a-z]?)(.*)", label)
        if match is None:
            raise ValueError(f"unsupported POTCAR dataset label: {label}")
        symbol, suffix = match.groups()
        if symbol in setups and setups[symbol] != suffix:
            raise ValueError(f"multiple POTCAR setups for {symbol}: {setups[symbol]!r} and {suffix!r}")
        setups[symbol] = suffix
    return setups


def _potcar_element(label: str) -> str:
    match = re.fullmatch(r"([A-Z][a-z]?)(.*)", label)
    if match is None:
        raise ValueError(f"unsupported POTCAR dataset label: {label}")
    return match.group(1)


def validate_vca_configuration(
    *,
    potcar: str | Path,
    vca_weights: Sequence[float] | None,
    virtual_symbol: str,
    components: Sequence[str],
) -> None:
    """Validate the one-site VASP-VCA contract before a run.

    VARNEB currently represents one physical site by coincident component atoms.
    The source POTCAR and its VCA vector must therefore contain exactly one
    contiguous component block; all other POTCAR datasets retain weight one.
    This rejects a mismatched POTCAR/INCAR pair before an expensive calculation.
    """

    component_list = [str(component) for component in components]
    if len(component_list) < 2:
        raise ValueError("VCA needs at least two coincident component datasets")
    if component_list[0] != virtual_symbol:
        raise ValueError(
            "the physical virtual symbol must be the first VCA component "
            f"({virtual_symbol!r} != {component_list[0]!r})"
        )
    if vca_weights is None:
        raise ValueError("VCA calculator requested but the source INCAR has no VCA tag")

    labels = potcar_dataset_labels(potcar)
    elements = [_potcar_element(label) for label in labels]
    weights = np.asarray(vca_weights, dtype=float).reshape(-1)
    if len(weights) != len(elements):
        raise ValueError(
            "VCA weight count must equal the number of POTCAR datasets "
            f"({len(weights)} != {len(elements)})"
        )
    if not np.isfinite(weights).all():
        raise ValueError("VCA weights must be finite")

    blocks = [
        start
        for start in range(len(elements) - len(component_list) + 1)
        if elements[start : start + len(component_list)] == component_list
    ]
    if len(blocks) != 1:
        raise ValueError(
            "POTCAR must contain exactly one contiguous VCA component block "
            f"for {component_list}; found {len(blocks)}"
        )
    start = blocks[0]
    component_indices = np.arange(start, start + len(component_list))
    component_weights = weights[component_indices]
    if (component_weights < 0.0).any() or not np.isclose(component_weights.sum(), 1.0, atol=1e-10):
        raise ValueError("VCA component weights must be nonnegative and sum to one")
    noncomponent = np.ones(len(weights), dtype=bool)
    noncomponent[component_indices] = False
    if not np.allclose(weights[noncomponent], 1.0, atol=1e-10):
        raise ValueError("all non-virtual POTCAR datasets must have VCA weight one")


class ExplicitPotcarVasp(Vasp):
    """ASE VASP calculator that restores the exact licensed source POTCAR."""

    def __init__(self, *, source_potcar: str | Path, **kwargs):
        self.source_potcar = Path(source_potcar).resolve()
        super().__init__(**kwargs)

    def write_input(self, atoms, properties=None, system_changes=None):
        super().write_input(atoms, properties=properties, system_changes=system_changes)
        shutil.copy2(self.source_potcar, Path(self.directory) / "POTCAR")


def expand_virtual_site(
    atoms: Atoms,
    *,
    virtual_symbol: str,
    components: Sequence[str],
) -> tuple[Atoms, list[list[int]]]:
    """Expand one physical virtual site into coincident VASP components."""

    indices = [index for index, symbol in enumerate(atoms.get_chemical_symbols()) if symbol == virtual_symbol]
    if len(indices) != 1:
        raise ValueError(f"expected exactly one {virtual_symbol} virtual site, found {len(indices)}")
    virtual_index = indices[0]
    symbols: list[str] = []
    scaled_positions = []
    physical_to_expanded: list[list[int]] = []
    original_scaled = atoms.get_scaled_positions(wrap=False)
    for index, (symbol, position) in enumerate(zip(atoms.get_chemical_symbols(), original_scaled)):
        mapped = []
        new_symbols = list(components) if index == virtual_index else [symbol]
        for new_symbol in new_symbols:
            mapped.append(len(symbols))
            symbols.append(new_symbol)
            scaled_positions.append(position)
        physical_to_expanded.append(mapped)
    expanded = Atoms(symbols, scaled_positions=scaled_positions, cell=atoms.cell, pbc=atoms.pbc)
    return expanded, physical_to_expanded


class VirtualCrystalCalculator(Calculator):
    """Expose a VASP VCA calculation as a physical, non-overlapping Atoms model."""

    implemented_properties = ["energy", "free_energy", "forces", "stress"]

    def __init__(
        self,
        base_calculator: Calculator,
        *,
        virtual_symbol: str,
        components: Sequence[str],
    ):
        super().__init__()
        self.base_calculator = base_calculator
        self.virtual_symbol = virtual_symbol
        self.components = tuple(components)
        self.directory = getattr(base_calculator, "directory", None)

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        if atoms is None:
            raise ValueError("VirtualCrystalCalculator requires atoms")
        expanded, mapping = expand_virtual_site(
            atoms,
            virtual_symbol=self.virtual_symbol,
            components=self.components,
        )
        expanded.calc = self.base_calculator
        energy = float(expanded.get_potential_energy())
        forces = np.asarray(expanded.get_forces(), dtype=float)
        stress = np.asarray(expanded.get_stress(), dtype=float)
        physical_forces = np.asarray([forces[indices].sum(axis=0) for indices in mapping])
        self.results = {"energy": energy, "forces": physical_forces, "stress": stress}
        free_energy = expanded.calc.results.get("free_energy")
        if free_energy is not None:
            self.results["free_energy"] = float(free_energy)


def vasp_input_fingerprints(source_dir: str | Path) -> dict:
    """Fingerprint the licensed VASP inputs copied to every VCNEB image."""

    source = Path(source_dir)
    records = {}
    for name in ("INCAR", "KPOINTS", "POTCAR"):
        path = source / name
        if not path.is_file():
            raise FileNotFoundError(path)
        contents = path.read_bytes()
        records[name] = {"path": str(path.resolve()), "bytes": len(contents), "sha256": hashlib.sha256(contents).hexdigest()}
    return records


def collect_vasp_params(calc: Vasp) -> dict:
    params = {}
    for name in [
        "float_params",
        "exp_params",
        "string_params",
        "int_params",
        "bool_params",
        "list_int_params",
        "list_bool_params",
        "list_float_params",
        "special_params",
        "dict_params",
        "input_params",
    ]:
        values = getattr(calc, name, None)
        if not values:
            continue
        for key, value in values.items():
            if value is not None:
                params[key] = value
    return params


def read_vasp_input_params(source_dir: str | Path) -> tuple[dict, dict, Path]:
    source = Path(source_dir)
    incar = source / "INCAR"
    kpoints = source / "KPOINTS"
    potcar = source / "POTCAR"
    for path in [incar, kpoints, potcar]:
        if not path.exists():
            raise FileNotFoundError(path)

    incar_reader = Vasp()
    incar_reader.read_incar(str(incar))
    incar_params = collect_vasp_params(incar_reader)
    for runtime_key in ("directory", "command", "txt", "label", "atoms"):
        incar_params.pop(runtime_key, None)

    kpoints_reader = Vasp()
    kpoints_reader.read_kpoints(str(kpoints))
    kpoint_params = {
        key: value
        for key, value in collect_vasp_params(kpoints_reader).items()
        if key in {"kpts", "gamma", "reciprocal", "kpts_nintersections"}
    }
    return incar_params, kpoint_params, potcar


def validate_vasp_static_parameters(parameters: Mapping) -> dict:
    """Validate the VASP image contract required by manager-owned VCNEB.

    A VASP image must return energy, forces and stress for the coordinates set
    by VARNEB.  ``relax``-like VASP settings would update atoms or the cell a
    second time and invalidate the NEB force evaluation, so they are rejected
    before an external executable is called.
    """

    normalized = {str(key).lower(): value for key, value in dict(parameters).items()}
    for key, expected in REQUIRED_VCNEB_STATIC_PARAMETERS.items():
        observed = normalized.get(key)
        try:
            matches = int(observed) == expected
        except (TypeError, ValueError):
            matches = False
        if not matches:
            raise ValueError(
                f"VASP VCNEB images require {key.upper()}={expected}, got {observed!r}"
            )
    return normalized


def prepare_vasp_static_parameters(
    source_dir: str | Path,
    *,
    overrides: Optional[Mapping] = None,
) -> tuple[dict, Path]:
    """Read an endpoint input and return the static VASP parameters for images."""

    incar_params, kpoint_params, potcar = read_vasp_input_params(source_dir)
    params = {str(key).lower(): value for key, value in incar_params.items()}
    params.update({str(key).lower(): value for key, value in kpoint_params.items()})
    params.update(
        {
            **REQUIRED_VCNEB_STATIC_PARAMETERS,
            "lcharg": False,
            "lwave": False,
        }
    )
    if overrides:
        params.update({str(key).lower(): value for key, value in dict(overrides).items()})
    try:
        params["setups"] = potcar_setups(potcar)
    except ValueError as exc:
        # ``--validate-only`` regression fixtures intentionally use a small
        # non-PAW placeholder.  A real VASP run will still reject that file.
        if "no TITEL records" not in str(exc):
            raise
    return validate_vasp_static_parameters(params), potcar


def attach_vasp_calculators(
    images: list[Atoms],
    *,
    source_dir: str | Path,
    workdir: str | Path,
    command: str,
    overrides: Optional[Mapping] = None,
    vca_virtual_symbol: Optional[str] = None,
    vca_components: Optional[Sequence[str]] = None,
) -> None:
    params, potcar = prepare_vasp_static_parameters(source_dir, overrides=overrides)

    root = Path(workdir)
    root.mkdir(parents=True, exist_ok=True)
    for image_index, image in enumerate(images):
        image_dir = root / f"{image_index:02d}"
        image_dir.mkdir(parents=True, exist_ok=True)
        write(image_dir / "POSCAR.start", image, format="vasp", direct=True, vasp5=True)
        shutil.copy2(potcar, image_dir / "POTCAR")
        base = ExplicitPotcarVasp(
            source_potcar=potcar,
            directory=str(image_dir),
            command=command,
            txt="vasp.out",
            **params,
        )
        if vca_virtual_symbol is None and vca_components is None:
            image.calc = base
        elif vca_virtual_symbol is not None and vca_components:
            validate_vca_configuration(
                potcar=potcar,
                vca_weights=params.get("vca"),
                virtual_symbol=vca_virtual_symbol,
                components=vca_components,
            )
            # Validate the physical structure at attachment time, not after the
            # first costly VASP call inside the optimizer.
            expand_virtual_site(
                image,
                virtual_symbol=vca_virtual_symbol,
                components=vca_components,
            )
            image.calc = VirtualCrystalCalculator(
                base,
                virtual_symbol=vca_virtual_symbol,
                components=vca_components,
            )
        else:
            raise ValueError("pass both vca_virtual_symbol and vca_components, or neither")


def default_vasp_command(ncores: int, executable: str) -> str:
    return os.environ.get("VASP_COMMAND", f"mpirun -np {ncores} {executable}")


def cached_vasp_static_endpoint_calculator(
    summary_path: str | Path,
    atoms: Atoms,
    *,
    endpoint: str,
    n_images: int,
    source_dir: str | Path,
    directory: str | Path,
) -> SinglePointCalculator:
    """Reuse a validated static endpoint so distributed workers remain interior-only."""

    if endpoint not in {"initial", "final"}:
        raise ValueError("endpoint must be 'initial' or 'final'")
    source = Path(summary_path).resolve()
    payload = json.loads(source.read_text(encoding="utf-8"))
    expected_index = 0 if endpoint == "initial" else n_images - 1
    expected_mode = f"fixed_{endpoint}_endpoint_static_scf"
    if payload.get("status") != "completed" or payload.get("execution_mode") != expected_mode:
        raise ValueError(f"{source} is not a completed {expected_mode} result")
    if payload.get("evaluated_image_index") != expected_index or payload.get("n_images") != n_images:
        raise ValueError(f"{source} does not match the requested {n_images}-image {endpoint} endpoint")
    recorded = ((payload.get("endpoint_structures") or {}).get(endpoint) or {}).get("sha256")
    if recorded != endpoint_structure_record(atoms).get("sha256"):
        raise ValueError(f"{source} structure does not match the requested {endpoint} endpoint")
    fingerprints = payload.get("licensed_input_fingerprints")
    if not isinstance(fingerprints, Mapping):
        raise ValueError(f"{source} lacks VASP input fingerprints")
    for name, active in vasp_input_fingerprints(source_dir).items():
        cached = fingerprints.get(name)
        if not isinstance(cached, Mapping) or cached.get("sha256") != active["sha256"]:
            raise ValueError(f"{source} {name} fingerprint does not match the active VASP input")
    forces = np.asarray(payload.get("forces_eV_per_A"), dtype=float)
    stress = np.asarray(payload.get("stress_eV_per_A3_voigt"), dtype=float)
    if forces.shape != (len(atoms), 3) or stress.shape != (6,):
        raise ValueError(f"{source} contains invalid force or stress arrays")
    calculator = SinglePointCalculator(atoms, energy=float(payload["potential_energy_eV"]), forces=forces, stress=stress)
    calculator.directory = str(Path(directory).resolve())
    return calculator
