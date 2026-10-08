"""Actual same-center G0 data, unit Jacobian and honest coverage limits."""

import copy
import json

import numpy as np
import pytest

from scripts.analyze_hfo2_joint_probe_reuse import EVIDENCE
from scripts.analyze_hfo2_restricted_quadratic import analyze_restricted, main


def evidence():
    manifest = json.loads((EVIDENCE / "work_probe_manifest.json").read_text())
    results = {i: json.loads((EVIDENCE / "point_audits/calculations" / f"{i:02d}" / "point_audit.json").read_text()) for i in range(8)}
    return manifest, results


def test_measured_coupling_and_nonstationary_release_coverage():
    report = analyze_restricted(*evidence())
    first = report["release_records"][0]
    assert first["softening_eV_A2"] == pytest.approx(2.3577754436039178)
    assert first["softening_fraction"] == pytest.approx(.1236356116081476)
    assert first["released_cell_curvature_eV_A2"] == pytest.approx(16.712583111960833)
    assert first["released_atomic_offset_A"] == pytest.approx(.021648157278814723)
    assert first["release_energy_change_at_zero_cell_coordinate_eV_cell"] == pytest.approx(-.0010945686515419764)
    assert report["sampling_coverage"]["all_released_coordinates_outside_probe_convex_hull"]
    assert all(not x["within_axis_probe_convex_hull"] for x in report["response_samples"])
    assert report["sampling_coverage"]["minimum_release_L1_A_over_interval"] == pytest.approx(first["released_atomic_offset_A"])
    assert not report["other_reference_Hessians_merged"]
    assert not report["unsampled_stability_known"]
    assert not report["conditional_DFT_surface_established"]
    assert not report["independent_barrier_prediction_validated"]
    assert not report["saddle_certified"]
    assert report["clamped_unmeasured_atomic_dimension"] == 35
    assert report["clamped_unmeasured_cell_dimension"] == 5
    assert report["new_DFT_calls"] == 0


def test_short_step_predicts_existing_long_axis_data_not_new_holdouts():
    report = analyze_restricted(*evidence())
    assert [x["index"] for x in report["retrospective_axis_crosschecks"]] == [2, 3, 6, 7]
    assert report["maximum_axis_energy_residual_eV_cell"] == pytest.approx(4.169314664268474e-5)
    assert report["maximum_axis_projected_gradient_residual_eV_A"] == pytest.approx(.004427275369460575)
    assert report["long_step_data_seen_before_model"]
    assert report["stability_screen_uses_two_step_data"]
    assert not report["floor_is_total_DFT_error_bound"]


def test_explicit_dimensionless_strain_jacobian_and_input_preservation():
    manifest, points = evidence()
    originals = copy.deepcopy((manifest, points))
    report = analyze_restricted(manifest, points)
    assert (manifest, points) == originals
    j = np.sqrt(2.) * manifest["cell_scale_A"]
    first = report["release_records"][0]
    assert report["cell_coordinate_per_biaxial_strain_A"] == pytest.approx(j)
    assert first["restricted_released_biaxial_strain_curvature_eV_cell"] == pytest.approx(first["released_cell_curvature_eV_A2"] * j * j)
    for point in report["response_samples"]:
        assert point["biaxial_strain"] * j == pytest.approx(point["cell_coordinate_A"])


@pytest.mark.parametrize("problem", ["wrong_cell_direction", "atomic_cell_mixture", "formula_units", "center_energy", "missing"])
def test_modified_coordinate_or_evidence_contract_is_refused(problem):
    manifest, results = evidence()
    if problem == "wrong_cell_direction":
        manifest["directions"][1][-6:] = [0., 0., 1., 0., 0., 0.]
    elif problem == "atomic_cell_mixture":
        manifest["directions"][0][-1] = .1
    elif problem == "formula_units":
        manifest["formula_units"] = 1
    elif problem == "center_energy":
        manifest["center_energy_eV_cell"] = float("nan")
    else:
        del results[3]
    with pytest.raises(ValueError):
        analyze_restricted(manifest, results)


def test_existing_report_refused_without_raw_read(tmp_path, monkeypatch):
    path = tmp_path / "result.json"
    path.write_text("original")
    monkeypatch.setattr("sys.argv", ["analysis", "--raw-root", str(tmp_path / "absent"), "--output", str(path)])
    with pytest.raises(FileExistsError):
        main()
    assert path.read_text() == "original"
