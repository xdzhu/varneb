import pytest

from scripts.prepare_gan_600eV_fft_restart_canary import FFT_TAGS, grid_locked_incar


def test_grid_lock_only_appends_original_pilot_mesh():
    source = b" ENCUT = 600.000000\n PSTRESS = 457.0\n NSW = 1\n"
    changed = grid_locked_incar(source)
    assert changed == source + b"".join(f" {tag}\n".encode() for tag in FFT_TAGS)


def test_grid_lock_rejects_existing_override():
    with pytest.raises(ValueError):
        grid_locked_incar(b" ENCUT = 600.000000\n PSTRESS = 457.0\n"
                          b" NSW = 1\n NGZ = 36\n")
