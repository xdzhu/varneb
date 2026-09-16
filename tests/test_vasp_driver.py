"""No-DFT gateway checks for the VASP VCNEB driver."""

from __future__ import annotations

import json
import runpy
import subprocess
import sys
from pathlib import Path

from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read, write
from vcneb import endpoint_structure_record, interpolate_vcneb


ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / "examples" / "run_vcneb_vasp.py"


def _endpoint(directory: Path, atoms: Atoms) -> None:
    directory.mkdir()
    write(directory / "CONTCAR", atoms, format="vasp")
    (directory / "INCAR").write_text("ENCUT = 500\nIBRION = 2\nNSW = 60\nISIF = 3\n", encoding="utf-8")
    (directory / "KPOINTS").write_text(
        "Automatic mesh\n0\nGamma\n1 1 1\n0 0 0\n",
        encoding="utf-8",
    )
    # This is deliberately not a real PAW dataset: --validate-only must not
    # parse or execute it, merely prove that every image gets an isolated file.
    (directory / "POTCAR").write_text("test POTCAR -- no DFT\n", encoding="utf-8")


def _module():
    return runpy.run_path(str(DRIVER))


def test_vasp_driver_validate_only_writes_static_7_image_preflight(tmp_path: Path, monkeypatch) -> None:
    initial = tmp_path / "initial"
    final = tmp_path / "final"
    _endpoint(initial, Atoms("Ba", cell=[4, 4, 4], pbc=True))
    _endpoint(
        final,
        Atoms("Ba", scaled_positions=[[0.1, 0.0, 0.0]], cell=[4.1, 4, 4], pbc=True),
    )
    workdir = tmp_path / "run"
    monkeypatch.setenv("VCNEB_GIT_REVISION", "remote-sync-test-vasp")
    result = subprocess.run(
        [
            sys.executable,
            str(DRIVER),
            "--initial",
            str(initial),
            "--final",
            str(final),
            "--workdir",
            str(workdir),
            "--validate-only",
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
    assert set(payload["licensed_input_fingerprints"]) == {"INCAR", "KPOINTS", "POTCAR"}
    assert all(len(item["sha256"]) == 64 for item in payload["licensed_input_fingerprints"].values())
    assert {key: payload["calculator_parameters"][key] for key in ("ibrion", "nsw", "isif", "isym")} == {
        "ibrion": -1,
        "nsw": 0,
        "isif": 2,
        "isym": 0,
    }
    assert len(payload["calculator_reports"]) == 7
    assert payload["git_revision"] == "remote-sync-test-vasp"
    assert all((workdir / f"{index:02d}" / "POTCAR").exists() for index in range(7))


def test_vasp_static_only_evaluates_fixed_initial_endpoint_once(tmp_path: Path, monkeypatch) -> None:
    module = _module()
    initial, final = tmp_path / "initial", tmp_path / "final"
    _endpoint(initial, Atoms("Ba", cell=[4, 4, 4], pbc=True))
    _endpoint(final, Atoms("Ba", scaled_positions=[[0.1, 0.0, 0.0]], cell=[4.1, 4, 4], pbc=True))
    calls: list[int] = []

    def fake_attach(images, *, source_dir, workdir, command, overrides):
        for index, atoms in enumerate(images):
            calls.append(index)
            calculator = SinglePointCalculator(
                atoms,
                energy=-8.0,
                forces=[[0.02, 0.0, 0.0]],
                stress=[0.1, 0.2, 0.3, 0.0, 0.0, 0.0],
            )
            calculator.directory = str(Path(workdir) / f"{index:02d}")
            atoms.calc = calculator

    monkeypatch.setattr(sys, "argv", [
        str(DRIVER), "--initial", str(initial), "--final", str(final),
        "--workdir", str(tmp_path / "static"), "--static-only",
    ])
    monkeypatch.setitem(module["main"].__globals__, "attach_vasp_calculators", fake_attach)
    module["main"]()

    summary = json.loads((tmp_path / "static" / "vasp_static_summary.json").read_text(encoding="utf-8"))
    assert calls == list(range(7))
    assert summary["execution_mode"] == "fixed_initial_endpoint_static_scf"
    assert summary["evaluated_image_index"] == 0
    assert summary["potential_energy_eV"] == -8.0
    assert summary["max_force_eV_per_A"] == 0.02


def test_vasp_static_only_can_evaluate_fixed_final_endpoint(tmp_path: Path, monkeypatch) -> None:
    module = _module()
    initial, final = tmp_path / "initial", tmp_path / "final"
    _endpoint(initial, Atoms("Ba", cell=[4, 4, 4], pbc=True))
    _endpoint(final, Atoms("Ba", scaled_positions=[[0.1, 0.0, 0.0]], cell=[4.1, 4, 4], pbc=True))

    def fake_attach(images, *, source_dir, workdir, command, overrides):
        for index, atoms in enumerate(images):
            calculator = SinglePointCalculator(atoms, energy=-7.0, forces=[[0.01, 0.0, 0.0]], stress=[0.0] * 6)
            calculator.directory = str(Path(workdir) / f"{index:02d}")
            atoms.calc = calculator

    monkeypatch.setattr(sys, "argv", [
        str(DRIVER), "--initial", str(initial), "--final", str(final),
        "--workdir", str(tmp_path / "static-final"), "--static-only", "--static-endpoint", "final",
    ])
    monkeypatch.setitem(module["main"].__globals__, "attach_vasp_calculators", fake_attach)
    module["main"]()
    summary = json.loads((tmp_path / "static-final" / "vasp_static_summary.json").read_text(encoding="utf-8"))
    assert summary["execution_mode"] == "fixed_final_endpoint_static_scf"
    assert summary["evaluated_image_index"] == 6


def test_vasp_final_static_records_the_mapped_path_endpoint_not_raw_input(tmp_path: Path, monkeypatch) -> None:
    module = _module()
    initial, final = tmp_path / "initial", tmp_path / "final"
    initial_atoms = Atoms("O2", scaled_positions=[[0.2, 0.0, 0.0], [0.7, 0.0, 0.0]], cell=[4, 4, 4], pbc=True)
    # The same physical sites in the opposite same-species order make the
    # distinction observable: automatic NEB mapping reorders the endpoint.
    final_atoms = Atoms("O2", scaled_positions=[[0.7, 0.0, 0.0], [0.2, 0.0, 0.0]], cell=[4, 4, 4], pbc=True)
    _endpoint(initial, initial_atoms)
    _endpoint(final, final_atoms)

    def fake_attach(images, *, source_dir, workdir, command, overrides):
        for index, atoms in enumerate(images):
            calculator = SinglePointCalculator(atoms, energy=-7.0, forces=[[0.01, 0.0, 0.0]] * len(atoms), stress=[0.0] * 6)
            calculator.directory = str(Path(workdir) / f"{index:02d}")
            atoms.calc = calculator

    monkeypatch.setattr(sys, "argv", [
        str(DRIVER), "--initial", str(initial), "--final", str(final),
        "--workdir", str(tmp_path / "static-final-mapped"), "--static-only", "--static-endpoint", "final",
        "--mapping", "auto",
    ])
    monkeypatch.setitem(module["main"].__globals__, "attach_vasp_calculators", fake_attach)
    module["main"]()

    summary = json.loads((tmp_path / "static-final-mapped" / "vasp_static_summary.json").read_text(encoding="utf-8"))
    expected = interpolate_vcneb(
        read(initial / "CONTCAR"),
        read(final / "CONTCAR"),
        n_images=7,
        align_cells=True,
        mic=True,
        cell_interpolation="log_strain",
        mapping="auto",
        align_translation=True,
        minimum_distance=1.6,
        maximum_deformation=0.10,
    )[-1]
    assert summary["endpoint_structures"]["final"]["sha256"] == endpoint_structure_record(expected)["sha256"]
    assert summary["endpoint_structures"]["final"]["sha256"] != endpoint_structure_record(read(final / "CONTCAR"))["sha256"]
