import copy
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("migration_gate", ROOT / "scripts/validate_vasp_migration_gate.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def probe():
    return {"mode": "converged_static_contract_probe", "all_passed": True, "isym": -1, "symprec": 1e-4,
            "records": [{"structure": "/case/gan_image15.vasp", "status": "converged_static_passed",
                         "repeated_input_write_passed": True, "energy_eV": -21.9, "maximum_force_eV_per_A": 0.13}]}


def test_migration_compare_converged_results_with_explicit_scope():
    old = probe(); new = copy.deepcopy(old)
    new["records"][0]["structure"] = "/new/gan_image15.vasp"
    new["records"][0]["energy_eV"] += 1e-6
    result = module.compare_probes(old, new)
    assert result["status"] == "ok"
    assert "not_full_array_or_path_certification" in result["scope"]
    new["records"][0]["energy_eV"] += 1e-2
    assert module.compare_probes(old, new)["status"] == "failed"


@pytest.mark.parametrize("change", ["initialization", "policy", "nonfinite", "write", "empty", "duplicate"])
def test_migration_rejects_invalid_or_incompatible_evidence(change):
    old = probe(); new = copy.deepcopy(old)
    if change == "initialization": new["mode"] = "initialization_probe_only_not_converged_scf"
    elif change == "policy": new["symprec"] = 1e-12
    elif change == "nonfinite": new["records"][0]["energy_eV"] = float("nan")
    elif change == "write": new["records"][0]["repeated_input_write_passed"] = False
    elif change == "empty": new["records"] = []
    elif change == "duplicate": new["records"] *= 2
    with pytest.raises(ValueError): module.compare_probes(old, new)


def test_hf_gan_layout_and_fixed_policy():
    text = (ROOT / "cluster/hf_gan_vcneb_vasp_distributed.slurm").read_text()
    assert "#SBATCH --ntasks=291" in text and "#SBATCH --ntasks-per-node=97" in text
    assert "--image-workers 9" in text and "--ntasks=32" in text
    assert "--exclusive --exact --nodes=1" in text and "--mpi=pmi2" in text
    assert "--no-climb" in text and "--resume-trajectory" in text
    assert "--vasp-symprec 1e-4" in text and "--fmax 0.10" in text
