"""Signed shape-error accounting; no clipping or invented material observations."""
from pathlib import Path
import importlib.util
import json
from copy import deepcopy
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "matched_frozen_release_tests",
    ROOT
    / "benchmarks/hfo2_channels/20261008/G3_frozen_release_E173_20261015/analyze.py",
)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
CASE = ROOT / "benchmarks/hfo2_channels/20261008/G3_frozen_release_E173_20261015"


@pytest.fixture(scope="module")
def material():
    summary = json.loads((CASE / "controls_summary.json").read_bytes())
    native = json.loads(
        (CASE / "admission_complete_native_preflight.json").read_bytes()
    )
    return summary, native, m.analyze(summary, native)


@pytest.mark.parametrize("seed", range(6))
def test_shape_error_decomposition_and_common_energy_offset_invariance(seed):
    rng = np.random.default_rng(seed)
    v = rng.normal(size=6)
    out = m.decompose_shape(
        *v[:4], model_released_origin=v[4], dft_released_origin=v[5]
    )
    assert out["released_shape_error_meV_cell"] == pytest.approx(
        out["frozen_shape_error_meV_cell"]
        + out["missed_release_shape_change_meV_cell"],
        abs=1e-9,
    )
    assert abs(out["error_decomposition_identity_defect_meV_cell"]) < 1e-9
    # All supplied quantities are already differences from one original center.
    # Adding the same absolute energy zero to original raw totals changes none.
    assert out["DFT_release_change_relative_origin_meV_cell"] == pytest.approx(
        ((v[3] - v[2]) - v[5]) * 1000
    )
    shifted = (v + 10000.0) - 10000.0
    again = m.decompose_shape(
        *shifted[:4], model_released_origin=shifted[4], dft_released_origin=shifted[5]
    )
    for name in out:
        if isinstance(out[name], float):
            assert again[name] == pytest.approx(out[name], abs=1e-8)


def test_higher_released_candidate_energy_is_preserved_not_clipped_or_called_invalid():
    out = m.decompose_shape(
        0.1, 0.08, 0.12, 0.15, model_released_origin=-0.01, dft_released_origin=-0.02
    )
    assert out["DFT_released_minus_frozen_meV_cell"] == pytest.approx(30.0)
    assert out["DFT_release_change_relative_origin_meV_cell"] == pytest.approx(50.0)
    assert out["signed_DFT_release_not_clipped"]


def test_missing_or_nonfinite_values_are_not_repaired():
    with pytest.raises(ValueError):
        m.decompose_shape(
            0.0, 0.0, np.nan, 0.0, model_released_origin=0.0, dft_released_origin=0.0
        )


@pytest.mark.parametrize(
    "fault",
    [
        "incomplete",
        "wrong_job",
        "unparsed",
        "source_mismatch",
        "control_not_required",
        "physical_holdout",
        "false_certificate",
    ],
)
def test_partial_or_incompatible_native_proof_is_rejected_before_any_material_loading(
    fault, monkeypatch
):
    summary = dict(
        status="eight_frozen_R_statics_completed_not_PES_or_independent_forecast",
        job_id="28865656",
        new_SCFs=8,
        total_attempts=8,
        rows=[dict(unit=i) for i in range(8)],
        source_original_centers_and_released_origins_reused=True,
        independent_material_prediction_or_released_stability_or_TS_certified=False,
        physical_reserved_geometry_generated_or_read=False,
    )
    native = dict(
        controls=dict(
            job_id="28865656",
            original_native_SCFs_reparsed=8,
            sources=deepcopy(summary["rows"]),
        ),
        original_runtime_files=380,
        original_native_SCFs=12,
        old_cache_hits_before_first_novel=7,
        controls_required=True,
        new_DFT=0,
    )
    if fault == "incomplete":
        summary["rows"].pop()
    elif fault == "wrong_job":
        summary["job_id"] = "another"
    elif fault == "unparsed":
        native["controls"]["original_native_SCFs_reparsed"] = 7
    elif fault == "source_mismatch":
        native["controls"]["sources"][0]["unit"] = 99
    elif fault == "control_not_required":
        native["controls_required"] = False
    elif fault == "physical_holdout":
        summary["physical_reserved_geometry_generated_or_read"] = True
    else:
        summary[
            "independent_material_prediction_or_released_stability_or_TS_certified"
        ] = True
    monkeypatch.setattr(
        m.faces,
        "load_inputs",
        lambda *a: pytest.fail("must reject incomplete proof before loading"),
    )
    with pytest.raises(ValueError, match="complete exact"):
        m.analyze(summary, native)


def test_all_eight_actual_controls_both_steps_and_signed_pointwise_identity_are_retained(
    material,
):
    summary, native, out = material
    assert (
        m.sha(CASE / "controls_summary.json")
        == "f49c662d705d4733dbcf1d78bb4e0260a0165df4c0a0df0af3e468f56fc8ffc8"
    )
    assert (
        m.sha(CASE / "admission_complete_native_preflight.json")
        == "8f0635eecad789677d0313c759924f2665a5f05186b334e9cb84ef8afde37d3f"
    )
    assert (
        native["controls"]["sources"] == summary["rows"] and len(summary["rows"]) == 8
    )
    assert [(r["channel"], r["curvature_probe_step_A"]) for r in out["rows"]] == [
        ("preserving", 0.01),
        ("preserving", 0.02),
        ("T", 0.01),
        ("T", 0.02),
    ]
    for r in out["rows"]:
        assert [p["point_id"] for p in r["points"]] == [1, 2, 3, 4]
        for p in r["points"]:
            assert p["released_shape_error_meV_cell"] == pytest.approx(
                p["frozen_shape_error_meV_cell"]
                + p["missed_release_shape_change_meV_cell"],
                abs=1e-8,
            )
    assert out["completed_static_attempts"] == 8
    assert out["completed_static_call_proxy_core_hours"] == pytest.approx(
        8.13681406982243, abs=1e-10
    )
    assert not out["physical_error_stability_TS_or_forecast_advantage_certified"]
    assert out["comparison_is_not_full_B2_B5_or_T_Cmma_control_matrix"]


@pytest.mark.parametrize(
    "i,ef,er,rel",
    [
        (0, 0.07588016876398755, 4.671403794771295, 4.606959080774686),
        (1, 0.11653548121566028, 4.72075845972527, 4.606959080774686),
        (2, 0.4162635209878788, 1.3271772093247154, 1.1687612932291813),
        (3, 0.4031595155926308, 1.3460706421835107, 1.1687612932291813),
    ],
)
def test_actual_matched_cohort_summary_statistics(material, i, ef, er, rel):
    r = material[-1]["rows"][i]
    assert r["max_abs_frozen_shape_error_meV_cell"] == pytest.approx(ef, abs=1e-9)
    assert r["max_abs_released_shape_error_meV_cell"] == pytest.approx(er, abs=1e-9)
    assert r["max_abs_DFT_release_shape_change_meV_cell"] == pytest.approx(
        rel, abs=1e-9
    )


@pytest.mark.parametrize(
    "fault", ["missing_point", "control", "input", "released_value"]
)
def test_changed_actual_control_cannot_be_silently_compared(material, fault):
    summary, native, _ = material
    summary = deepcopy(summary)
    native = deepcopy(native)
    if fault == "missing_point":
        summary["rows"][0]["point_id"] = 99
    elif fault == "control":
        summary["rows"][0]["coordinates_full_A"][-1] = 1e-6
    elif fault == "input":
        summary["rows"][0]["native_validation"]["input_sha256"]["INPUT"] = "0" * 64
    else:
        summary["rows"][0]["matched_released_energy_eV_cell"] += 0.001
    native["controls"]["sources"] = deepcopy(summary["rows"])
    with pytest.raises(ValueError):
        m.analyze(summary, native)


def test_actual_launch_receipts_bind_complete_controls_and_at_most_two_allocations():
    receipts = json.loads((CASE / "matched_arm_submission_receipts.json").read_bytes())
    proof = json.loads((CASE / "launch_receipt.json").read_bytes())
    assert receipts["accepted"]["job_id"] == proof["matched_arm_job"] == "28866132"
    assert proof["old_static_controls_sacct"].split("|")[1:3] == ["COMPLETED", "0:0"]
    assert receipts["submission_intent"]["complete_native_preflight_sha256"] == m.sha(
        CASE / "admission_complete_native_preflight.json"
    )
    assert len(receipts["held_validation"]["admission"]["committed_G3_queue"]) == 2
    assert receipts["submission_intent"]["admission"]["control"] == dict(
        job_id="28865656", state_exit=["COMPLETED", "0:0"]
    )
    assert proof["additional_training_reservation"] == 0


def test_recorded_parent_alias_is_accepted_only_for_the_same_real_native_directory(
    tmp_path,
):
    spec = importlib.util.spec_from_file_location(
        "readonly_first_native_path_test", CASE / "audit_first.py"
    )
    audit = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(audit)
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    evaluation = tmp_path / "arm/eval-000000-unit"
    (evaluation / "scf_000000").mkdir(parents=True)
    raw = bundle / "../arm/eval-000000-unit/scf_000000"
    directory = audit.recorded_evaluation_directory(
        str(raw), evaluation, tmp_path / "arm"
    )
    assert str(directory / "scf_000000") == str(raw)
    assert directory.resolve() == evaluation.resolve()
    with pytest.raises(ValueError):
        audit.recorded_evaluation_directory(str(raw), evaluation, tmp_path / "bundle")
    other = tmp_path / "arm/eval-other"
    (other / "scf_000000").mkdir(parents=True)
    with pytest.raises(ValueError):
        audit.recorded_evaluation_directory(str(raw), other, tmp_path / "arm")
