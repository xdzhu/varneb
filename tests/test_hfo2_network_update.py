"""Real frozen material identities and figure integrity, not new DFT tests."""

import copy
import json
from pathlib import Path

import numpy as np
import pytest
from ase.io import read

pytest.importorskip("pymatgen")

from scripts.analyze_hfo2_network_update import analyze_case, structure_audit
import scripts.plot_hfo2_network_update as plotter


ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / "benchmarks/hfo2_channels/20261008/network_update_20261009"


@pytest.fixture(scope="module")
def actual_case():
    return analyze_case(CASE / "network_specification.json")


def test_actual_phase_sweep_and_nonpolar_centres_are_not_stationary_claims(actual_case):
    assert actual_case["existing_image_records"] == 37
    assert actual_case["new_DFT_calls"] == 0
    assert not actual_case["TS_certified"] and not actual_case["full_G1_passed"]
    assert [c["ordinary_residual_passed"] for c in actual_case["channels"]] == [True, False, False, False]
    preserving, reversing = actual_case["channels"][-2:]
    for candidate, expected in ((preserving, "Pbcn"), (reversing, "Pbca")):
        peak = candidate["rows"][candidate["peak_view_image_index"]]
        assert [s["symbol"] for s in peak["structure_audit"]["symmetry_sweep"]] == [expected] * 3
        assert [s["symprec_A"] for s in peak["structure_audit"]["symmetry_sweep"]] == [.001, .01, .05]
        assert not peak["structure_audit"]["standardized_wrapped_reordered_or_relaxed"]
        assert peak["NEB_cell_max_vector_eV_A"] > .10
    assert preserving["peak_energy_relative_PO_meV_fu"] == pytest.approx(65.27856176671776)
    assert reversing["peak_energy_relative_PO_meV_fu"] == pytest.approx(410.7436409649381)
    # Zero triplet coverage is demonstrably not a cubic-phase certificate.
    rpeak = reversing["rows"][reversing["peak_view_image_index"]]
    assert rpeak["triplet_captured_squared_norm_fraction"] < 1e-20
    assert abs(rpeak["Q_pattern_x_A"]) < 1e-8
    assert all(r["energy_relative_PO_meV_fu"] < actual_case["channels"][0]["rows"][-1]["energy_relative_PO_meV_fu"] for r in preserving["rows"])


def test_same_initial_and_thermodynamic_reverse_do_not_recompute_force(actual_case):
    t = actual_case["channels"][0]
    assert t["source_direction"] == "reverse"
    assert t["rows"][0]["source_image_index"] == 9
    assert t["rows"][-1]["source_image_index"] == 0
    assert t["rows"][0]["s_normalized"] == 0.
    assert t["rows"][-1]["s_normalized"] == 1.
    for c in actual_case["channels"]:
        assert c["rows"][0]["energy_relative_PO_meV_fu"] == pytest.approx(0., abs=1e-6)
        assert c["rows"][0]["NEB_max_vector_eV_A"] is None
        assert c["rows"][-1]["NEB_max_vector_eV_A"] is None
        assert np.all(np.diff([r["s_normalized"] for r in c["rows"]]) > 0)


def test_structure_audit_preserves_coordinates_and_reports_elemental_occupancy():
    atom = read(ROOT / "benchmarks/hfo2_channels/20261008/reference_variants/T.vasp")
    cell, frac = atom.cell.array.copy(), atom.get_scaled_positions(wrap=False).copy()
    result = structure_audit(atom)
    np.testing.assert_array_equal(atom.cell.array, cell)
    np.testing.assert_array_equal(atom.get_scaled_positions(wrap=False), frac)
    assert result["ordered"] and not result["oxidation_states_guessed"]
    assert [s["symbol"] for s in result["symmetry_sweep"]] == ["P4_2/nmc"] * 3


@pytest.mark.parametrize("problem", ["not_periodic", "wrong_order", "left_handed", "nan", "collision"])
def test_invalid_structure_rejected_before_interpreting_symmetry(problem):
    atom = read(ROOT / "benchmarks/hfo2_channels/20261008/reference_variants/T.vasp")
    if problem == "not_periodic":
        atom.pbc = False
    elif problem == "wrong_order":
        atom = atom[[4, 1, 2, 3, 0, 5, 6, 7, 8, 9, 10, 11]]
    elif problem == "left_handed":
        cell = atom.cell.array.copy(); cell[0] *= -1
        atom.set_cell(cell)
    elif problem == "nan":
        atom.positions[0, 0] = np.nan
    else:
        atom.positions[1] = atom.positions[0]
    with pytest.raises(ValueError):
        structure_audit(atom)


def test_frozen_material_report_is_replayed_before_drawing(actual_case):
    data = plotter.build_data(CASE / "network_specification.json", CASE / "material_observations_v3.json")
    assert len(data["rows"]) == 37
    assert data["all_four_energy_profiles_in_source_data"]
    assert not data["reversing_energy_profile_in_panel_a"]
    assert not data["TS_or_final_channel_ranking_certified"]
    fig = plotter.make_figure(data)
    assert len(fig.axes) == 4 and len(fig.legends) == 1
    positions = [a.get_position().bounds for a in fig.axes]
    for i, a in enumerate(fig.axes):
        assert not a.get_title() and all(s.get_visible() for s in a.spines.values())
        label = next(t for t in a.texts if t.get_text() == f"({chr(97+i)})")
        assert label.get_fontweight() == "normal"
        assert label.get_fontsize() == 13
    assert positions[0][1] == positions[1][1] and positions[2][1] == positions[3][1]
    assert positions[0][0] == positions[2][0] and positions[1][0] == positions[3][0]
    legend_box = fig.legends[0].get_window_extent(fig.canvas.get_renderer())
    assert all(legend_box.y0 > a.get_window_extent().y1 for a in fig.axes)
    import matplotlib.pyplot as plt
    plt.close(fig)


@pytest.mark.parametrize("fault", ["phase", "energy", "force", "hash"])
def test_stale_material_fields_fail_closed_before_figure(tmp_path, monkeypatch, actual_case, fault):
    changed = copy.deepcopy(actual_case)
    row = changed["channels"][2]["rows"][4]
    if fault == "phase":
        row["structure_audit"]["symmetry_sweep"][0]["symbol"] = "Pm-3m"
    elif fault == "energy":
        row["energy_relative_PO_meV_fu"] += 1.
    elif fault == "force":
        row["NEB_max_vector_eV_A"] += .01
    else:
        changed["T_reference_sha256"] = "0" * 64
    path = tmp_path / "stale.json"
    path.write_text(json.dumps(changed))
    monkeypatch.setattr(plotter, "analyze_case", lambda _: actual_case)
    with pytest.raises(ValueError):
        plotter.build_data(CASE / "network_specification.json", path)


def test_existing_output_refused_before_material_read(tmp_path, monkeypatch):
    path = tmp_path / "keep"; path.mkdir()
    monkeypatch.setattr("sys.argv", ["plot", "--specification", "missing", "--material-report", "missing", "--output", str(path)])
    with pytest.raises(FileExistsError):
        plotter.main()


def test_warning_cache_metadata_is_preserved_without_weakening_scientific_fields(actual_case):
    recorded = copy.deepcopy(actual_case)
    warnings = recorded["channels"][2]["rows"][4]["structure_audit"]["warnings"]
    warnings.append("recorded cache-dependent warning")
    before = copy.deepcopy(recorded)
    differences = plotter.verify_material_replay(recorded, actual_case)
    assert recorded == before
    assert len(differences) == 1
    assert differences[0]["field"] == "/channels/2/rows/4/structure_audit/warnings"
    assert differences[0]["recorded"][-1] == "recorded cache-dependent warning"
    recorded["channels"][2]["rows"][4]["structure_audit"]["symmetry_sweep"][0]["number"] = 225
    with pytest.raises(ValueError):
        plotter.verify_material_replay(recorded, actual_case)
