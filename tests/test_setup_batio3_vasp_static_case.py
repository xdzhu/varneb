"""Regression checks for licensed-POTCAR BTO VASP case preparation."""

from __future__ import annotations

import json
import runpy
from pathlib import Path

from ase import Atoms
from ase.io import write

from vcneb.provenance import endpoint_structure_record


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "setup_batio3_vasp_static_case.py"


def _reference(path: Path, initial: Atoms, final: Atoms) -> None:
    path.write_text(json.dumps({"endpoint_structures": {"initial": endpoint_structure_record(initial), "final": endpoint_structure_record(final)}}), encoding="utf-8")


def test_setup_bto_vasp_static_case_copies_only_approved_inputs_after_identity_check(tmp_path: Path) -> None:
    module = runpy.run_path(str(SCRIPT))
    initial = Atoms("BaTiO3", cell=[4, 4, 4], pbc=True)
    final = Atoms("BaTiO3", scaled_positions=[[0.0, 0.0, 0.0]] * 5, cell=[4.1, 4.1, 4.1], pbc=True)
    initial_path, final_path = tmp_path / "initial.vasp", tmp_path / "final.vasp"
    write(initial_path, initial, format="vasp")
    write(final_path, final, format="vasp")
    reference = tmp_path / "reference.json"
    _reference(reference, initial, final)
    potcar, incar, kpoints = tmp_path / "POTCAR", tmp_path / "INCAR", tmp_path / "KPOINTS"
    potcar.write_text("licensed POTCAR bytes\n", encoding="utf-8")
    incar.write_text("ENCUT = 600\n", encoding="utf-8")
    kpoints.write_text("Gamma\n", encoding="utf-8")
    destination = tmp_path / "case"

    result = module["setup"](
        initial_path=initial_path,
        final_path=final_path,
        reference_identity=reference,
        potcar=potcar,
        destination=destination,
        incar=incar,
        kpoints=kpoints,
    )
    assert result["status"] == "prepared_no_dft"
    assert (destination / "initial" / "POTCAR").read_bytes() == potcar.read_bytes()
    assert (destination / "final" / "CONTCAR").exists()
    assert not (destination / "final" / "POTCAR").exists()
    assert json.loads((destination / "setup_manifest.json").read_text(encoding="utf-8"))["endpoint_comparison"]["matches"] is True


def test_setup_bto_vasp_static_case_rejects_nonidentical_endpoint(tmp_path: Path) -> None:
    module = runpy.run_path(str(SCRIPT))
    initial = Atoms("BaTiO3", cell=[4, 4, 4], pbc=True)
    final = Atoms("BaTiO3", cell=[4.1, 4.1, 4.1], pbc=True)
    initial_path, final_path = tmp_path / "initial.vasp", tmp_path / "final.vasp"
    write(initial_path, initial, format="vasp")
    write(final_path, final, format="vasp")
    reference = tmp_path / "reference.json"
    _reference(reference, initial, initial)
    for name in ("POTCAR", "INCAR", "KPOINTS"):
        (tmp_path / name).write_text(name, encoding="utf-8")

    try:
        module["setup"](
            initial_path=initial_path,
            final_path=final_path,
            reference_identity=reference,
            potcar=tmp_path / "POTCAR",
            destination=tmp_path / "case",
            incar=tmp_path / "INCAR",
            kpoints=tmp_path / "KPOINTS",
        )
    except ValueError as exc:
        assert "do not match" in str(exc)
    else:
        raise AssertionError("mismatched endpoint was accepted")
