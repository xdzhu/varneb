import pytest

from scripts.prepare_gan_600eV_basin_uninterrupted import long_incar


def test_long_run_keeps_physics_and_enables_wave_checkpoint():
    source = (b" ENCUT = 600.000000\n PSTRESS = 457.0\n ISIF = 3\n"
              b" IBRION = 2\n NSW = 10\n LWAVE = .FALSE.\n")
    result = long_incar(source)
    assert result == (source.replace(b" NSW = 10\n", b" NSW = 100\n")
                            .replace(b" LWAVE = .FALSE.\n", b" LWAVE = .TRUE.\n"))


def test_long_run_rejects_different_encut():
    with pytest.raises(ValueError):
        long_incar(b" ENCUT = 1000.000000\n PSTRESS = 457.0\n ISIF = 3\n"
                   b" IBRION = 2\n NSW = 10\n LWAVE = .FALSE.\n")
