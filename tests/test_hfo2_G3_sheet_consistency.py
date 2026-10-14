"""Actual training arithmetic and explicit synthetic controls, no new DFT."""
from copy import deepcopy
import importlib.util
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / "benchmarks/hfo2_channels/20261008/G3_sheet_consistency_E172_20261015"
spec = importlib.util.spec_from_file_location(
    "sheet_diagnostic_tests", CASE / "analyze.py"
)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


@pytest.fixture(scope="module")
def sheet():
    data, rows, report = m.load_T()
    return data, rows, report, m.diagnose(data, rows)


def test_all_original_nodes_edges_loops_and_raw_labels_are_retained(sheet):
    _, rows, report, result = sheet
    assert len(rows) == 25 and report["new_SCFs_reparsed"] == 118
    assert (
        len(result["nodes"]) == 25
        and len(result["edges"]) == 40
        and len(result["loops"]) == 16
    )
    s = result["summary"]
    assert s["relative_to_released_origin_energy_range_meV_cell"] == pytest.approx(
        [-38.530500440174364, 0.0], abs=1e-10
    )
    assert s["total_radius_range_A"] == pytest.approx(
        [0.04579226130057099, 0.4637805580466008], abs=1e-12
    )
    assert s["minimum_original_radius_margin_A"] == pytest.approx(
        0.03621944195339921, abs=1e-12
    )
    assert s["released_gradient_range_eV_A"] == pytest.approx(
        [0.00652512162716755, 0.029784325041568457], abs=1e-12
    )
    assert s["new_DFT_calls"] == 0
    for r in result["edges"]:
        assert r["prescribed_control_endpoint_work_meV_cell"] == 0.0
        assert (
            abs(r["released_endpoint_work_meV_cell"])
            <= r["endpoint_released_work_Cauchy_bound_meV_cell"] + 2e-9
        )


def test_retained_and_full_work_diagnostics_are_distinct_not_error_bars(sheet):
    r = sheet[-1]
    s = r["summary"]
    ret = s["largest_retained_only_energy_work_difference"]
    assert ret["ids"] == [14, 7]
    assert ret["retained_only_energy_minus_work_meV_cell"] == pytest.approx(
        1.2091245416555734, abs=1e-9
    )
    assert ret["released_endpoint_work_meV_cell"] == pytest.approx(
        1.0938736182258546, abs=1e-9
    )
    assert ret["full_energy_minus_trapezoid_work_meV_cell"] == pytest.approx(
        0.11525092342971541, abs=1e-9
    )
    assert s["largest_full_endpoint_quadrature_difference"]["ids"] == [3, 9]
    assert s["largest_full_endpoint_quadrature_difference"][
        "full_energy_minus_trapezoid_work_meV_cell"
    ] == pytest.approx(1.0394056615714766, abs=1e-9)
    for edge in r["edges"]:
        assert edge["retained_only_energy_minus_work_meV_cell"] == pytest.approx(
            edge["released_endpoint_work_meV_cell"]
            + edge["full_energy_minus_trapezoid_work_meV_cell"],
            abs=1e-9,
        )


def test_closed_energy_telescope_is_not_zero_endpoint_gradient_circulation(sheet):
    loops = sheet[-1]["loops"]
    assert all(r["raw_closed_loop_energy_change_meV_cell"] == 0.0 for r in loops)
    r = sheet[-1]["summary"]["largest_retained_loop_circulation"]
    assert r["ids"] == [7, 3, 9, 16, 7]
    assert r["retained_endpoint_circulation_meV_cell"] == pytest.approx(
        -2.000383803111482, abs=1e-9
    )
    assert r["released_endpoint_circulation_meV_cell"] == pytest.approx(
        1.4627714851978684, abs=1e-9
    )
    assert r["full_endpoint_circulation_meV_cell"] == pytest.approx(
        -0.5376123179136171, abs=1e-9
    )
    for r in loops:
        assert r["retained_endpoint_circulation_meV_cell"] + r[
            "released_endpoint_circulation_meV_cell"
        ] == pytest.approx(r["full_endpoint_circulation_meV_cell"], abs=1e-9)


def toy_node(x, g, energy, identifier):
    x = np.asarray(x)
    g = np.asarray(g)
    return dict(
        id=identifier,
        q=x[:2],
        x=x,
        full=np.r_[x, np.zeros(4)],
        gf=np.r_[g, np.zeros(4)],
        g=g,
        energy=energy,
        residual=float(np.linalg.norm(g[2:])),
    )


def test_quadratic_synthetic_control_exact_energy_work_with_nonzero_released_work():
    a = np.eye(36)[:, :2]
    r = np.eye(36)[:, 2:]
    h = np.diag(np.arange(1, 37, dtype=float))
    h[0, 2] = h[2, 0] = 0.3
    x = np.zeros(36)
    y = x.copy()
    x[:3] = [0.1, 0.2, 0.25]
    y[:3] = [0.2, 0.1, -0.15]
    left = toy_node(x, h @ x, 0.5 * x @ h @ x, 0)
    right = toy_node(y, h @ y, 0.5 * y @ h @ y, 1)
    out = m.segment(left, right, a, r)
    assert out["full_energy_minus_trapezoid_work_meV_cell"] == pytest.approx(
        0.0, abs=1e-10
    )
    assert abs(out["released_endpoint_work_meV_cell"]) > 1.0
    assert out["retained_only_energy_minus_work_meV_cell"] == pytest.approx(
        out["released_endpoint_work_meV_cell"], abs=1e-10
    )
    rr = np.linalg.qr(np.random.default_rng(4).normal(size=(34, 34)))[0]
    other = m.segment(left, right, a, r @ rr)
    for k in out:
        if isinstance(out[k], float):
            assert other[k] == pytest.approx(out[k], abs=1e-10)


def test_zero_endpoint_residual_does_not_bound_unsampled_energy_error():
    # E(r)=3r²-2r³ has E(1)-E(0)=1, but both endpoint derivatives are zero.
    x = np.zeros(36)
    y = x.copy()
    y[2] = 1.0
    out = m.segment(
        toy_node(x, np.zeros(36), 0.0, 0),
        toy_node(y, np.zeros(36), 1.0, 1),
        np.eye(36)[:, :2],
        np.eye(36)[:, 2:],
    )
    assert out["released_endpoint_work_meV_cell"] == 0.0
    assert out["endpoint_released_work_Cauchy_bound_meV_cell"] == 0.0
    assert out["full_energy_minus_trapezoid_work_meV_cell"] == 1000.0
    assert out["endpoint_Cauchy_bound_is_not_an_energy_error_bound"]


@pytest.mark.parametrize(
    "fault",
    ["missing", "extra", "label", "control", "gradient", "nonfinite", "embedding"],
)
def test_invalid_or_unregistered_candidate_records_are_refused(sheet, fault):
    data, rows, _, _ = sheet
    rows = deepcopy(rows)
    if fault == "missing":
        rows.pop()
    elif fault == "extra":
        rows.append(rows[-1])
    elif fault == "label":
        rows[1]["q_A"] = [0.3, 0.3]
    elif fault == "control":
        rows[1]["coordinates_full_A"][-1] = 1e-6
    elif fault == "gradient":
        rows[1]["released_gradient_norm_eV_A"] = 0.0
    elif fault == "nonfinite":
        rows[1]["energy_eV_cell"] = float("nan")
    else:
        rows[1]["coordinates_internal_A"][0] += 0.01
    with pytest.raises(ValueError):
        m.diagnose(data, rows)


def test_unmatched_native_report_bytes_are_rejected_before_material_loading(
    tmp_path, monkeypatch
):
    p = tmp_path / "fake.json"
    p.write_text("{}")
    monkeypatch.setattr(m, "T_AUDIT", p)
    monkeypatch.setattr(
        m.faces, "validate", lambda *a: pytest.fail("must reject first")
    )
    with pytest.raises(ValueError, match="exact original"):
        m.load_T()


def test_complete_analysis_owns_no_physical_accuracy_or_new_data_claim():
    out = m.analyze()
    assert not out["source_native_SCFs_newly_reparsed"]
    assert not out[
        "physical_error_bound_or_stability_or_unique_branch_or_TS_or_independent_prediction_certified"
    ]
    assert not out["physical_reserved_geometry_generated_or_read"]
    assert out["additional_training_points"] == out["scheduler_writes"] == 0
    assert (
        out["formula_units"] == 4
        and out["source_pressure_GPa"] == out["source_physical_strain"] == 0.0
    )
    assert len(out["source_six_physical_inputs"]) == 6
