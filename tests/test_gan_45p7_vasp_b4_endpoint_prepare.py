"""B4 endpoint refinement must preserve the 600-eV electronic contract."""

import shutil

import pytest

from scripts import prepare_gan_45p7_vasp_b4_endpoint_refine as preparer


def _source(tmp_path, monkeypatch):
    source = tmp_path / "original_static"
    source.mkdir()
    archived = preparer.Path(
        "paper/VARNEB_CPC/evidence/gan_vasp_endpoint_statics_20260930"
    )
    for filename in ("INCAR", "KPOINTS", "POSCAR", "OUTCAR"):
        shutil.copy2(archived / f"initial_00_{filename}", source / filename)
    (source / "POTCAR").write_bytes(b"licensed-potential-test-placeholder")
    monkeypatch.setitem(preparer.SOURCE_HASHES, "POTCAR",
                        preparer.sha256(source / "POTCAR"))
    return source


def test_b4_refinement_changes_only_ionic_cell_settings(tmp_path, monkeypatch):
    source = _source(tmp_path, monkeypatch)
    original_incar = (source / "INCAR").read_bytes()
    work = tmp_path / "isolated_work"
    manifest = preparer.prepare(source, work)
    generated = (work / "case" / "INCAR").read_bytes()
    assert b" ENCUT = 600.000000\n" in generated
    assert b" SYMPREC = 1.00e-04\n" in generated
    assert b" ISYM = -1\n" in generated
    assert b" IBRION = 2\n" in generated
    assert b" ISIF = 3\n" in generated
    assert b" PSTRESS = 457.0\n" in generated
    assert b" EDIFFG = -0.02\n" in generated
    assert (source / "INCAR").read_bytes() == original_incar
    assert preparer.sha256(source / "KPOINTS") == preparer.sha256(
        work / "case" / "KPOINTS")
    assert preparer.sha256(source / "POTCAR") == preparer.sha256(
        work / "case" / "POTCAR")
    assert manifest["status"] == "inputs_finalized_no_DFT"
    with pytest.raises(FileExistsError):
        preparer.prepare(source, work)


def test_b4_refinement_rejects_changed_original_endpoint(tmp_path, monkeypatch):
    source = _source(tmp_path, monkeypatch)
    (source / "POSCAR").write_text("changed", encoding="utf-8")
    work = tmp_path / "isolated_work"
    with pytest.raises(ValueError, match="original B4 source changed"):
        preparer.prepare(source, work)
    assert not work.exists()
