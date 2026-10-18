"""Completed native branch vs connected polygon, without a new DFT call."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT/'benchmarks/hfo2_channels/20261008/preserving_pattern_branches_E272_20261018'
spec = importlib.util.spec_from_file_location('P_native_connection_test', CASE/'analyze_prefix.py')
reader = importlib.util.module_from_spec(spec); spec.loader.exec_module(reader)


@pytest.fixture(scope='module')
def analysis():
    return reader.analyze()


def test_original_native_nine_image_chain_is_reproduced(analysis):
    report, _ = analysis; original = report['original']
    assert original['replayed_fmax_eV_A'] == pytest.approx(.0988168085874911, abs=1e-12)
    assert original['sampled_polygon_peak_meV_fu'] == pytest.approx(96.82818300279905, abs=1e-9)
    assert original['ordinary_residual_pass'] and original['interior_local_peak_images'] == [4]


def test_lowered_center_changes_neighbors_tangents_and_does_not_pass_NEB(analysis):
    report, _ = analysis; original = report['original']; new = report['replaced_central_image']
    assert new['sampled_polygon_peak_meV_fu'] == pytest.approx(85.49988851655144, abs=1e-9)
    assert new['sampled_peak_image'] == 3 and new['interior_local_peak_images'] == [3, 5]
    assert not new['ordinary_residual_pass']
    assert new['replayed_fmax_eV_A'] == pytest.approx(.1256114037176295, abs=1e-12)
    assert new['per_interior_image_fmax_eV_A'][3] < .02
    assert new['per_interior_image_fmax_eV_A'][2] > .12
    assert new['per_interior_image_fmax_eV_A'][4] > .12
    for index in (0, 1, 5, 6):
        assert new['per_interior_image_fmax_eV_A'][index] == original['per_interior_image_fmax_eV_A'][index]
    assert new['geometry_gate']['adjacent_segment_cosines'][3] < 0


def test_physical_replay_preserves_matched_controls_and_separate_secants(analysis):
    report, _ = analysis; branch = report['signed_branch']
    assert branch['local_native_SCFs_reparsed'] == 4
    assert branch['arms']['A']['dE_minus_F_meV_fu'] == pytest.approx(-13.59925840177256, abs=1e-8)
    assert branch['arms']['J']['dE_minus_F_meV_fu'] == pytest.approx(-13.617081922348, abs=1e-8)
    for row in branch['arms'].values():
        assert row['released_gradient_norm_eV_A'] < .03 and row['retained_pattern_defect_A'] < 1e-10
        assert not row['full_orbital_pseudopotential_bytes_rechecked_locally']
    secants = branch['directional_secants_at_A_eV_A2']
    assert secants['gradient'] > 1 and .9 < secants['energy'] < 1
    assert not branch['full_stability_at_A_or_J_certified']


def test_only_center_replaced_endpoints_order_and_continuous_lift_preserved(analysis):
    report, images = analysis; originals, _ = reader.original_images()
    assert len(images) == 9 and report['n_free_interior_images'] == 7
    for index in range(9):
        assert images[index].get_chemical_symbols() == ['Hf']*4+['O']*8
        if index != 4:
            np.testing.assert_array_equal(images[index].positions, originals[index].positions)
            np.testing.assert_array_equal(images[index].cell.array, originals[index].cell.array)
    assert np.any(np.asarray(report['fixed_per_atom_integer_lift']) != 0)
    assert report['new_DFT_calls'] == 0 and report['original_labels_and_forecast_scores_unchanged']
    assert not report['whole_E272_job_completed_here']
    assert not report['stationary_TS_or_MEP_or_corrected_barrier_certified']
    assert not report['interpolating_segments_evaluated']


def test_lift_rejects_a_nearby_but_different_frozen_structure():
    chart, _, _, _ = reader.setup(); images, _ = reader.original_images()
    wrong = chart.reference_atoms.copy(); wrong.positions[0, 0] += .01
    with pytest.raises(ValueError, match='differs from the original ordered peak'):
        reader.fixed_lift(images[4], wrong, chart.reference_atoms)


def test_lift_is_transport_not_new_minimum_image_search():
    chart, _, _, _ = reader.setup(); images, _ = reader.original_images()
    candidate = chart.reference_atoms.copy()
    candidate.positions[0, 0] += .001
    candidate.calc = reader.SinglePointCalculator(candidate, energy=0., forces=np.zeros((12, 3)), stress=np.zeros(6))
    lifted, integers = reader.fixed_lift(images[4], chart.reference_atoms, candidate)
    np.testing.assert_allclose(lifted.get_scaled_positions(wrap=False),
                               candidate.get_scaled_positions(wrap=False)+integers, atol=1e-14, rtol=0)
    assert reader.same_ordered_geometry(lifted, candidate)


def test_changed_archive_cannot_be_accepted(tmp_path, monkeypatch):
    import shutil
    shutil.copy(CASE/'prefix_0_minus_receipt.json', tmp_path/'prefix_0_minus_receipt.json')
    (tmp_path/'prefix_0_minus.tar.gz').write_bytes(b'not the native archive')
    monkeypatch.setattr(reader, 'CASE', tmp_path)
    with pytest.raises(ValueError, match='archive changed'):
        reader.restore(tmp_path/'out')
