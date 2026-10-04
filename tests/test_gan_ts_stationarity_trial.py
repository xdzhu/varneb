"""Guard the bounded same-contract GaN stationarity proposal."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from scripts.prepare_gan_600eV_ts_stationarity_trial import hybrid_newton_step


ROOT = Path(__file__).resolve().parents[1]
HESSIAN = ROOT / "paper/VARNEB_CPC/evidence/gan_45p7_split_20261004/joint_hessian.npz"
REPORT = ROOT / "benchmarks/numerical_integrity/gan_600eV_normal_strain_step_dependence_20261001.json"
TRIAL = ROOT / "benchmarks/stationarity/gan_45p7_20261004"


def test_correction_is_small_and_uses_measured_normal_energy_slopes():
    record = hybrid_newton_step(HESSIAN, REPORT)
    assert 0 < record["newton_shift_norm_A"] < 0.004
    assert np.count_nonzero(np.asarray(record["newton_shift_A"])) > 0
    for axis in (12, 13, 14):
        assert record["hybrid_energy_informed_gradient_eV_per_A"][axis] == (
            record["normal_energy_secants_eV_per_A"][str(axis)]
        )
        assert record["hybrid_energy_informed_gradient_eV_per_A"][axis] < -0.02


def test_wrong_hessian_bytes_are_rejected(tmp_path):
    wrong = tmp_path / "wrong.npz"
    wrong.write_bytes(b"not the audited local mode basis")
    with pytest.raises(ValueError, match="Hessian basis changed"):
        hybrid_newton_step(wrong, REPORT)


def test_stationarity_trial_outcome_is_hash_bound_and_not_overclaimed():
    manifest_bytes = (TRIAL / "manifest.json").read_bytes()
    manifest = json.loads(manifest_bytes)
    result = json.loads((TRIAL / "stationarity_audit.json").read_text())
    assert result["manifest_sha256"] == hashlib.sha256(manifest_bytes).hexdigest()
    assert manifest["electronic_contract"] == (
        "VASP6.3.2/PBE/Ga_d+N/600eV/Gamma8x8x6/EDIFF1e-7/ISYM-1/SYMPREC1e-4"
    )
    assert len(result["cases"]) == manifest["n_cases"] == 9
    assert not result["all_predeclared_trial_gates_pass"]
    assert result["center_max_atomic_force_eV_per_A"] < 0.01
    assert result["max_normal_energy_derived_pressure_residual_kbar"] < 2
    assert result["center_stress_residual_max_kbar"] > 2
    assert result["max_normal_energy_stress_mismatch_kbar"] > 2
    assert result["unstable_energy_curvature_eV_per_A2"] < 0
