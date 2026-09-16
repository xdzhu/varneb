"""Small VASP helpers for VC-NEB examples."""

from __future__ import annotations

import os
import hashlib
import json
import shutil
from pathlib import Path
from typing import Mapping, Optional

import numpy as np
from ase import Atoms
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
    return validate_vasp_static_parameters(params), potcar


def attach_vasp_calculators(
    images: list[Atoms],
    *,
    source_dir: str | Path,
    workdir: str | Path,
    command: str,
    overrides: Optional[Mapping] = None,
) -> None:
    params, potcar = prepare_vasp_static_parameters(source_dir, overrides=overrides)

    root = Path(workdir)
    root.mkdir(parents=True, exist_ok=True)
    for image_index, image in enumerate(images):
        image_dir = root / f"{image_index:02d}"
        image_dir.mkdir(parents=True, exist_ok=True)
        write(image_dir / "POSCAR.start", image, format="vasp", direct=True, vasp5=True)
        shutil.copy2(potcar, image_dir / "POTCAR")
        image.calc = Vasp(directory=str(image_dir), command=command, txt="vasp.out", **params)


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
