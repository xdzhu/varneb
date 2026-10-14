"""Completed0/+1% matrix, explicit sampled trends and unknown error bounds."""
import json

import numpy as np
import pytest

from scripts import audit_hfo2_clamped_network as network


@pytest.fixture(scope='module')
def result():
    return network.audit()


def test_all_eight_native_bands_share_exact_initial_cache_per_condition(result):
    assert result['ordinary_audited_training_paths'] == 8
    assert result['portable_raw_frame_records_reparsed'] == 72
    assert sorted(r['occurrences'] for r in result['initial_cache_groups']) == [4, 4]
    assert {r['substrate_strain'] for r in result['rows']} == {0., .01}
    for row in result['rows']:
        assert row['ordinary_fmax_eV_A'] <= .10
        assert row['forward_meV_fu']-row['reverse_meV_fu'] == pytest.approx(row['final_minus_initial_meV_fu'], abs=1e-10)
        assert row['barrier_error_bound_meV_fu'] is None
        assert len(row['raw_relative_energy_meV_fu']) == 9
    assert result['pressure_GPa'] == 0 and result['new_DFT_calls'] == result['scheduler_mutations'] == 0
    assert not result['new_mode_or_reference_selection']
    assert not result['unpublished_author_arrays_required']


def test_opposite_flip_responses_do_not_allow_cherrypicking_the_high_barrier(result):
    r = result['response']
    assert r['candidate_changes_meV_fu'] == pytest.approx(dict(M=-1.2304034162298194,
        T=12.347026513452874, preserving=28.5962481552815, reversing=-44.02433096583991), abs=1e-9)
    assert all(x['sampled_switching_candidate'] == 'preserving' for x in r['sampled_minima'].values())
    assert all(x['sampled_escape_candidate'] == 'T' for x in r['sampled_minima'].values())
    assert r['sampled_S_change_meV_fu'] == pytest.approx(-16.249221641828626, abs=1e-9)
    assert not r['raw_H1_favorable_signature'] and not r['H1_error_certified_rejection']
    assert not r['barrier_error_bounds_measured']
    assert not r['G3_modes_or_holdout_edges_selected']
    assert not r['physical_holdout_generated_or_read']
    assert all(not x['physical_ordering_resolved'] for x in r['sampled_minima'].values())


def test_all_original_reference_reports_include_four_complete_Cmma_frames():
    base = network.BASE/'clamped_plus1_reversing_terminal_E148_20261014/analysis'
    d = json.loads((base/'displacement.json').read_text())
    g = json.loads((base/'gradient.json').read_text())
    w = json.loads((base/'path_work.json').read_text())
    assert d['highest_image_index'] == g['sampled_peak_image'] == w['sampled_peak_image'] == 4
    assert d['T_triplet_rank3_fraction'][4] == pytest.approx(.16141026412529272)
    spans = g['sampled_peak']['representations']
    assert spans['T_geometric_triplet']['retained_squared_gradient_fraction'] < 1e-25
    assert spans['T_two_lowest_complete_optical_doublets']['retained_squared_gradient_fraction'] == pytest.approx(.2248412479727011)
    for frame in ('103', '106', '114', '127'):
        x = spans['Cmma_frame_'+frame]
        assert x['source_basis_rank'] == 4 and not x['source_modes_or_forces_ASR_corrected']
        assert x['retained_squared_gradient_fraction'] == pytest.approx(.63521999256, abs=1e-10)
    t = w['representations']['T_geometric_triplet']
    assert 250*t['cumulative_energy_minus_work'][4] == pytest.approx(35.85699318496349)
    assert abs(250*t['endpoint_energy_minus_work']) < .000002
    assert not w['G3_selection_or_holdout_access']


def groups():
    return {0.: {k: dict(forward_meV_fu=v) for k,v in [('M',4.),('T',3.),('preserving',2.),('reversing',5.)]},
            .01: {k: dict(forward_meV_fu=v) for k,v in [('M',5.),('T',4.),('preserving',1.),('reversing',6.)]}}


def test_a_favorable_sampled_signature_still_does_not_certify_H1():
    r = network.summarize(groups())
    assert r['raw_H1_favorable_signature']
    assert r['H1_primary_status'].startswith('not_established')
    assert not r['H1_error_certified_rejection'] and not r['barrier_error_bounds_measured']


@pytest.mark.parametrize('mode', ['condition', 'edge', 'negative', 'nonfinite', 'boolean'])
def test_incomplete_or_nonphysical_matrix_refused(mode):
    x = groups()
    if mode == 'condition': del x[.01]
    elif mode == 'edge': del x[0.]['reversing']
    else: x[0.]['M']['forward_meV_fu'] = {'negative':-1., 'nonfinite':np.nan, 'boolean':True}[mode]
    with pytest.raises(ValueError): network.summarize(x)
