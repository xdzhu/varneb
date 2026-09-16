"""No-DFT gateway checks for the VASP VCNEB driver."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from ase import Atoms
from ase.io import write


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
