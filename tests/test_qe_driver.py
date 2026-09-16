"""Static gateway checks for the QE VCNEB driver."""

from __future__ import annotations

import json
import runpy
import subprocess
import sys
from pathlib import Path

from ase import Atoms
from ase.io import write
import pytest


ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / "examples" / "run_vcneb_qe.py"


def _module():
    return runpy.run_path(str(DRIVER))


def test_qe_driver_defaults_to_7_total_images_and_0p10(monkeypatch) -> None:
    module = _module()
    monkeypatch.setattr(sys, "argv", [str(DRIVER), "--initial", "a", "--final", "b"])
    args = module["parse_args"]()
    assert args.n_images == 7
    assert args.fmax == 0.10
    assert args.kpts == [4, 4, 4]


def test_qe_driver_parameters_require_complete_species_and_static_cutoffs(monkeypatch) -> None:
    module = _module()
    monkeypatch.setattr(sys, "argv", [str(DRIVER), "--initial", "a", "--final", "b", "--pp", "Ba=Ba.upf"])
    args = module["parse_args"]()
    assert module["build_qe_parameters"](args, {"Ba"})["input_data"]["system"]["ecutwfc"] == 100.0
    with pytest.raises(ValueError, match="missing --pp"):
        module["build_qe_parameters"](args, {"Ba", "Ti"})
    args.ecutrho = 80.0
    with pytest.raises(ValueError, match="ecutrho"):
        module["build_qe_parameters"](args, {"Ba"})


def test_qe_driver_validate_only_writes_a_7_image_preflight(tmp_path) -> None:
    initial = Atoms("Ba", cell=[4, 4, 4], pbc=True)
    final = Atoms("Ba", scaled_positions=[[0.1, 0.0, 0.0]], cell=[4.1, 4, 4], pbc=True)
    initial_path, final_path = tmp_path / "initial.vasp", tmp_path / "final.vasp"
    write(initial_path, initial, format="vasp")
    write(final_path, final, format="vasp")
    (tmp_path / "Ba.upf").write_text('<UPF element="Ba" functional="PBE">\n', encoding="utf-8")
    workdir = tmp_path / "run"
    result = subprocess.run(
        [
            sys.executable, str(DRIVER), "--initial", str(initial_path), "--final", str(final_path),
            "--workdir", str(workdir), "--command", "pw.x", "--pseudo-dir", str(tmp_path),
            "--pp", "Ba=Ba.upf", "--validate-only",
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    assert "preflight passed" in result.stdout
    payload = json.loads((workdir / "vcneb_preflight.json").read_text(encoding="utf-8"))
    assert payload["n_images"] == 7
    assert payload["n_interior_images"] == 5
    assert payload["fmax_target_eV_per_A"] == 0.10
    assert payload["calculator_parameters"]["input_data"]["control"] == {
        "calculation": "scf", "tstress": True, "tprnfor": True
    }
    assert payload["pseudopotential_reports"][0]["sha256"]
