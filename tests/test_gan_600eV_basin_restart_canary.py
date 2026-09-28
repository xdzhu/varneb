from pathlib import Path

import pytest

from scripts.prepare_gan_600eV_basin_restart_canary import canary_incar, prepare


def test_canary_changes_only_nsw():
    raw = (b" ENCUT = 600.000000\n PSTRESS = 457.0\n ISIF = 3\n"
           b" IBRION = 2\n NSW = 10\n EDIFFG = -0.02\n POTIM = 0.25\n")
    changed = canary_incar(raw)
    assert changed == raw.replace(b" NSW = 10\n", b" NSW = 1\n")


@pytest.mark.parametrize("raw", [
    b" ENCUT = 1000.000000\n PSTRESS = 457.0\n ISIF = 3\n IBRION = 2\n NSW = 10\n",
    b" ENCUT = 600.000000\r\n PSTRESS = 457.0\r\n ISIF = 3\r\n IBRION = 2\r\n NSW = 10\r\n",
    b" ENCUT = 600.000000\n PSTRESS = 0.0\n ISIF = 3\n IBRION = 2\n NSW = 10\n",
])
def test_canary_rejects_drift(raw):
    with pytest.raises(ValueError):
        canary_incar(raw)


def test_prepare_requires_existing_audited_pilot(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        prepare(tmp_path / "absent", tmp_path / "absent.json", tmp_path / "result")
