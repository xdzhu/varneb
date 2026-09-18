"""Regression tests for VASP's static VCNEB image contract."""

from __future__ import annotations

import json
from pathlib import Path

from ase import Atoms
import pytest

from vcneb.provenance import endpoint_structure_record
from vcneb.vasp import (
    REQUIRED_VCNEB_STATIC_PARAMETERS,
    cached_vasp_static_endpoint_calculator,
    vasp_input_fingerprints,
    validate_vasp_static_parameters,
)


def test_vasp_static_contract_accepts_required_parameters_case_insensitively() -> None:
    parameters = {key.upper(): value for key, value in REQUIRED_VCNEB_STATIC_PARAMETERS.items()}
    validated = validate_vasp_static_parameters(parameters)
    assert validated == REQUIRED_VCNEB_STATIC_PARAMETERS


def test_vasp_static_contract_accepts_fully_disabled_symmetry_for_low_symmetry_images() -> None:
    parameters = dict(REQUIRED_VCNEB_STATIC_PARAMETERS, isym=-1)
    assert validate_vasp_static_parameters(parameters)["isym"] == -1


@pytest.mark.parametrize("key, value", [("ibrion", 2), ("nsw", 5), ("isif", 3), ("isym", 2)])
def test_vasp_static_contract_rejects_internal_relaxation(key, value) -> None:
    parameters = dict(REQUIRED_VCNEB_STATIC_PARAMETERS)
    parameters[key] = value
    with pytest.raises(ValueError, match=key.upper()):
        validate_vasp_static_parameters(parameters)


def test_cached_vasp_static_endpoint_reuses_only_matching_fingerprinted_input(tmp_path: Path) -> None:
    for name, text in {"INCAR": "ENCUT = 600\n", "KPOINTS": "Gamma\n", "POTCAR": "licensed bytes\n"}.items():
        (tmp_path / name).write_text(text, encoding="utf-8")
    atoms = Atoms("Ba", cell=[4, 4, 4], pbc=True)
    summary = {
        "status": "completed",
        "execution_mode": "fixed_initial_endpoint_static_scf",
        "evaluated_image_index": 0,
        "n_images": 7,
        "endpoint_structures": {"initial": endpoint_structure_record(atoms)},
        "licensed_input_fingerprints": vasp_input_fingerprints(tmp_path),
        "potential_energy_eV": -1.0,
        "forces_eV_per_A": [[0.1, 0.0, 0.0]],
        "stress_eV_per_A3_voigt": [0.0] * 6,
    }
    path = tmp_path / "summary.json"
    path.write_text(json.dumps(summary), encoding="utf-8")
    calculator = cached_vasp_static_endpoint_calculator(
        path, atoms, endpoint="initial", n_images=7, source_dir=tmp_path, directory=tmp_path / "00",
    )
    assert calculator.get_potential_energy(atoms) == -1.0
    (tmp_path / "POTCAR").write_text("changed bytes\n", encoding="utf-8")
    with pytest.raises(ValueError, match="POTCAR fingerprint"):
        cached_vasp_static_endpoint_calculator(
            path, atoms, endpoint="initial", n_images=7, source_dir=tmp_path, directory=tmp_path / "00",
        )
