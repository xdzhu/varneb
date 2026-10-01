"""A stricter B4 continuation cannot silently change the DFT contract."""

import shutil

import pytest

from scripts import prepare_gan_45p7_vasp_b4_endpoint_strict as strict


def _first_run(tmp_path, monkeypatch):
    archived = strict.Path("paper/VARNEB_CPC/evidence/gan_vasp_b4_endpoint_refine_20261001")
    source = tmp_path / "first_run"
    case = source / "case"
    case.mkdir(parents=True)
    shutil.copy2(archived / "manifest.json", source / "manifest.json")
    for name in ("INCAR", "KPOINTS", "POSCAR", "OUTCAR", "CONTCAR"):
        shutil.copy2(archived / name, case / name)
    (case / "POTCAR").write_bytes(b"licensed-potential-test-placeholder")
    monkeypatch.setitem(strict.EXPECTED_SOURCE_SHA256, "case/POTCAR",
                        strict.sha256(case / "POTCAR"))
    return source


def test_strict_continuation_changes_only_stop_settings(tmp_path, monkeypatch):
    source = _first_run(tmp_path, monkeypatch)
    work = tmp_path / "strict_run"
    manifest = strict.prepare(source, work)
    before = (source / "case/INCAR").read_bytes()
    after = (work / "case/INCAR").read_bytes()
    assert after == before.replace(b" NSW = 100\n", b" NSW = 30\n").replace(
        b" EDIFFG = -0.02\n", b" EDIFFG = -0.005\n"
    )
    assert manifest["electronic_contract"]["ENCUT_eV"] == 600
    for name in ("KPOINTS", "POTCAR"):
        assert strict.sha256(source / "case" / name) == strict.sha256(
            work / "case" / name)
    assert strict.sha256(source / "case/CONTCAR") == strict.sha256(
        work / "case/POSCAR")
    with pytest.raises(FileExistsError):
        strict.prepare(source, work)


def test_strict_continuation_rejects_first_run_drift(tmp_path, monkeypatch):
    source = _first_run(tmp_path, monkeypatch)
    (source / "case/OUTCAR").write_bytes(b"changed")
    work = tmp_path / "strict_run"
    with pytest.raises(ValueError, match="source changed"):
        strict.prepare(source, work)
    assert not work.exists()
