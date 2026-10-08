import copy
import json
from pathlib import Path

import numpy as np
import pytest

from scripts.analyze_hfo2_joint_probe_reuse import analyze, main


ROOT = Path(__file__).resolve().parents[1] / "benchmarks/hfo2_channels/20261008"


def evidence():
    manifest = json.loads((ROOT / "work_probe_manifest.json").read_text())
    results = {i: json.loads((ROOT / "point_audits/calculations" / f"{i:02d}" / "point_audit.json").read_text()) for i in range(8)}
    return manifest, results


def test_actual_probe_columns_retain_mixed_and_transverse_response():
    manifest, results = evidence()
    report = analyze(manifest, results)
    first, second = report["step_records"]
    np.testing.assert_allclose(first["raw_projected_eV_A2"],
                               [[4.671228720095983, -3.318399815220913],
                                [-3.3189817069456913, 19.070358555564752]], atol=1e-12)
    np.testing.assert_allclose(first["transverse_column_norms_eV_A2"],
                               [5.867122503121306, 7.094254351721581], atol=1e-12)
    assert first["projected_eigenvalues_eV_A2"][0] > 3.9
    assert second["reciprocity_relative_defect"] > first["reciprocity_relative_defect"]
    assert report["two_step_full_action_operator_spread_eV_A2"] == pytest.approx(.11787188791706943)
    assert report["center_full_gradient_norm_eV_A"] > .11
    assert not report["raw_source_audit_performed"]
    assert not report["full_joint_Hessian_measured"]
    assert not report["saddle_certified"]
    assert not report["conditional_surface_established"]
    assert not report["independent_barrier_prediction_validated"]
    historical = json.loads((ROOT / "work_probe_summary.json").read_text())
    for old in historical["pairs"]:
        current = report["step_records"][int(old["step_A"] == .02)]
        diagonal = current["raw_projected_eV_A2"][old["direction"]][old["direction"]]
        assert diagonal == pytest.approx(old["gradient_curvature_eV_A2"], abs=1e-12)


@pytest.mark.parametrize("change", ["missing", "duplicate_tag", "pressure", "bad_gradient", "bad_status", "bad_index"])
def test_incomplete_or_mismatched_evidence_is_refused(change):
    manifest, results = evidence()
    if change == "missing":
        del results[7]
    elif change == "duplicate_tag":
        manifest["points"][1] = copy.deepcopy(manifest["points"][0])
    elif change == "pressure":
        manifest["pressure_eV_A3"] = 1
    elif change == "bad_gradient":
        results[3]["gradient_eV_A"][0] = float("nan")
    elif change == "bad_status":
        results[1]["status"] = "pending"
    else:
        results[2]["index"] = 3
    with pytest.raises(ValueError):
        analyze(manifest, results)


def test_existing_report_is_refused_before_any_raw_source_read(tmp_path, monkeypatch):
    output = tmp_path / "existing.json"
    output.write_text("keep original receipt")
    monkeypatch.setattr("sys.argv", ["analyze", "--raw-root", str(tmp_path / "absent"), "--output", str(output)])
    with pytest.raises(FileExistsError):
        main()
    assert output.read_text() == "keep original receipt"
