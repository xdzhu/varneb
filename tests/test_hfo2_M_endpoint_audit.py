import json
from pathlib import Path

from ase.io import read
import numpy as np
import pytest

import scripts.audit_hfo2_M_endpoint as audit


def test_endpoint_reaudits_raw_force_and_stress_not_only_wrapper(tmp_path, monkeypatch):
    root = tmp_path / "M"
    source = root / "calculator/image_0000/scf_000000"
    source.mkdir(parents=True)
    (source / "INPUT").write_text("fixed fixture")
    (source / "STRU").write_text("mock geometry")
    (source / "OUT.ABACUS").mkdir()
    (source / "OUT.ABACUS/running_scf.log").write_text("synthetic only")
    monkeypatch.setattr(audit, "CONTRACT", {"INPUT": audit.sha256(source / "INPUT")})
    summary = {"status": "completed", "converged": True, "external_pressure_gpa": 0,
               "optimizer_steps": 0, "potential_energy_eV": -1}
    (root / "endpoint_relax_summary.json").write_text(json.dumps(summary))
    atoms = read(Path(__file__).resolve().parents[1] / "benchmarks/hfo2_channels/20261008/reference_variants/M_seed.vasp")
    monkeypatch.setattr(audit, "read", lambda *a, **k: atoms.copy())
    monkeypatch.setattr(audit, "symmetry_report", lambda a: {"symmetry": [{"number": 14}]})
    raw = {"energy": -1, "forces": np.zeros((12, 3)), "stress": np.zeros(6)}
    monkeypatch.setattr(audit, "audited_results", lambda p: raw)
    assert audit.audit(root, tmp_path / "pass.json")["n_real_audited_SCFs"] == 1
    raw["stress"][0] = 3/1602.176634
    with pytest.raises(ValueError, match="gates"):
        audit.audit(root, tmp_path / "fail.json")
