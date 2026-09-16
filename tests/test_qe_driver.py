"""Static gateway checks for the QE VCNEB driver."""

from __future__ import annotations

import json
import runpy
import subprocess
import sys
from pathlib import Path

from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
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


def test_qe_driver_rejects_mixing_approved_manifest_with_manual_pp(tmp_path) -> None:
    initial = Atoms("Ba", cell=[4, 4, 4], pbc=True)
    endpoint = tmp_path / "endpoint.vasp"
    write(endpoint, initial, format="vasp")
    manifest = tmp_path / "manifest.json"
    manifest.write_text('{"approval_status":"approved","species":{"Ba":{"filename":"Ba.UPF","md5":"00000000000000000000000000000000"}}}', encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(DRIVER), "--initial", str(endpoint), "--final", str(endpoint), "--command", "pw.x", "--pseudo-dir", str(tmp_path), "--pp-manifest", str(manifest), "--pp", "Ba=Ba.UPF", "--validate-only"],
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode != 0
    assert "either --pp-manifest" in result.stderr


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
    assert payload["endpoint_structures"]["initial"]["sha256"]
    assert payload["endpoint_structures"]["initial"]["n_atoms"] == 1
    assert payload["endpoint_structures"]["initial"]["sha256"] != payload["endpoint_structures"]["final"]["sha256"]
    assert payload["calculator_parameters"]["input_data"]["control"] == {
        "calculation": "scf", "tstress": True, "tprnfor": True
    }
    assert payload["pseudopotential_reports"][0]["sha256"]


def test_qe_driver_static_only_evaluates_fixed_initial_endpoint_once(tmp_path, monkeypatch) -> None:
    module = _module()
    initial = Atoms("Ba", cell=[4, 4, 4], pbc=True)
    final = Atoms("Ba", scaled_positions=[[0.1, 0.0, 0.0]], cell=[4.1, 4, 4], pbc=True)
    initial_path, final_path = tmp_path / "initial.vasp", tmp_path / "final.vasp"
    write(initial_path, initial, format="vasp")
    write(final_path, final, format="vasp")
    (tmp_path / "Ba.upf").write_text('<UPF element="Ba" functional="PBE">\\n', encoding="utf-8")

    calls: list[int] = []

    def fake_factory(*, parameters, command, pseudo_dir):
        def make(index, atoms, directory):
            calls.append(index)
            calculator = SinglePointCalculator(
                atoms,
                energy=-12.5,
                forces=[[0.01, 0.0, 0.0]],
                stress=[0.1, 0.2, 0.3, 0.0, 0.0, 0.0],
            )
            calculator.directory = str(directory)
            return calculator
        return make

    monkeypatch.setattr(sys, "argv", [
        str(DRIVER), "--initial", str(initial_path), "--final", str(final_path),
        "--workdir", str(tmp_path / "static"), "--command", "pw.x", "--pseudo-dir", str(tmp_path),
        "--pp", "Ba=Ba.upf", "--static-only",
    ])
    # runpy returns a result dictionary distinct from function.__globals__.
    # Patch the latter so this test can exercise the real driver flow without
    # invoking an external pw.x executable.
    monkeypatch.setitem(module["main"].__globals__, "make_ase_espresso_factory", fake_factory)
    module["main"]()

    summary = json.loads((tmp_path / "static" / "qe_static_summary.json").read_text(encoding="utf-8"))
    assert calls == list(range(7))  # immutable calculator directories remain auditable.
    assert summary["execution_mode"] == "fixed_initial_endpoint_static_scf"
    assert summary["evaluated_image_index"] == 0
    assert summary["potential_energy_eV"] == -12.5
    assert summary["max_force_eV_per_A"] == pytest.approx(0.01)
