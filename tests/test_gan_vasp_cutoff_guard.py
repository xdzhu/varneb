"""GaN benchmark must reject VASP cutoff drift before launching DFT."""

import pytest

from examples.run_vcneb_vasp import validate_expected_encut


def test_fixed_gan_cutoff_accepts_only_600_ev():
    validate_expected_encut({"encut": 600}, 600)
    with pytest.raises(ValueError, match="differs from the fixed benchmark"):
        validate_expected_encut({"encut": 1000}, 600)
    with pytest.raises(ValueError, match="differs from the fixed benchmark"):
        validate_expected_encut({"encut": 800}, 600)


def test_generic_vasp_run_can_use_its_own_explicit_contract():
    validate_expected_encut({"encut": 500}, None)
    with pytest.raises(ValueError, match="no valid ENCUT"):
        validate_expected_encut({}, 600)
    with pytest.raises(ValueError, match="finite and positive"):
        validate_expected_encut({"encut": 600}, -1)


def test_gan_launcher_requires_original_cutoff():
    from pathlib import Path

    launcher = Path(__file__).resolve().parents[1] / "cluster" / "hf_gan_vcneb_vasp_distributed.slurm"
    assert "--expected-encut-ev 600" in launcher.read_text(encoding="utf-8")
