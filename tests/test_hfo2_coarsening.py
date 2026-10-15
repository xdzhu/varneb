"""Elementary retrospective comparators, not a new physical predictor."""
from pathlib import Path
import importlib.util,itertools,json
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[1]
CASE=ROOT/'benchmarks/hfo2_channels/20261008/G3_coarsening_E184_20261015'
s=importlib.util.spec_from_file_location('E184_coarsening',CASE/'analyze.py')
a=importlib.util.module_from_spec(s);s.loader.exec_module(a)


@pytest.mark.parametrize('kind',['bilinear','tensor_quadratic'])
def test_constant_and_corner_interpolation_reproduces_values_and_unit_sum(kind):
    labels={ij:2. for ij in itertools.product((-2,0,2),repeat=2)}
    for ij in itertools.product(range(-2,3),repeat=2):
        p=a.interpolate(ij,labels,kind)
        assert p['value']==pytest.approx(2.,abs=1e-14)
        assert sum(v['weight'] for v in p['required_weights'])==pytest.approx(1.,abs=1e-14)


def test_tensor_quadratic_exactly_reproduces_mixed_fourth_degree_not_just_harmonic():
    def energy(i,j):return 1+i+2*j+i*j+3*i*i+2*j*j+4*i*i*j-2*i*j*j+5*i*i*j*j
    labels={ij:energy(*ij) for ij in itertools.product((-2,0,2),repeat=2)}
    for ij in itertools.product(range(-2,3),repeat=2):
        assert a.interpolate(ij,labels,'tensor_quadratic')['value']==pytest.approx(energy(*ij),abs=1e-12)


@pytest.mark.parametrize('kind',['bilinear','tensor_quadratic'])
def test_missing_corner_refuses_off_axis_but_zero_weight_does_not_invent_data(kind):
    labels={ij:float(ij[0]**2+ij[1]**2) for ij in itertools.product((-2,0,2),repeat=2) if 0 in ij}
    assert a.interpolate((1,1),labels,kind)['value'] is None
    assert a.interpolate((1,1),labels,kind)['missing_required_nodes']
    assert a.interpolate((1,0),labels,kind)['value'] is not None
    assert all(v['weight']!=0 for v in a.interpolate((1,0),labels,kind)['required_weights'])


def test_full_quadratic_fit_and_rank_deficient_plane_are_not_pseudoinverted():
    labels={ij:float(a.features(ij)@np.arange(6)) for ij in itertools.product((-2,0,2),repeat=2)}
    beta,rank=a.fit_quadratic(labels);assert rank==6
    np.testing.assert_allclose(beta,np.arange(6),atol=2e-12,rtol=0)
    beta,rank=a.fit_quadratic({ij:v for ij,v in labels.items() if 0 in ij})
    assert beta is None and rank==5


@pytest.mark.parametrize('index',[-3,3,.5])
def test_no_extrapolation_outside_registered_integer_grid(index):
    with pytest.raises(ValueError):a.one_axis(index,'tensor_quadratic')


@pytest.fixture(scope='module')
def analysis():return a.analyze()


def test_actual_split_uses_all_coarse_indices_and_keeps_holes_not_trials(analysis):
    t,p=analysis['rows']
    assert len(t['coarse_ids'])==len(p['coarse_ids'])==9
    assert len(t['available_coarse_ids'])==9 and len(t['observed_fine_ids'])==16
    assert len(p['available_coarse_ids'])==5 and len(p['observed_fine_ids'])==12
    assert p['quadratic_design_rank']==5 and p['quadratic_coefficients'] is None
    assert len(p['missing_coarse_ids'])==4 and len(p['original_fine_holes'])==4
    assert p['common_all_models_ids']==[]
    assert p['common_without_rank_deficient_fit_ids']==[1,2,3,4] and p['secondary_cohort_is_axis_only']
    for r in (t,p):
        assert r['coarse_geometry_and_branch_history_were_not_blinded_to_fine_preparation']
        assert r['coarse_nodes_started_from_fine_parents']
        assert set(r['observed_fine_ids']).isdisjoint(r['coarse_ids'])


def test_actual_source_cost_includes_calls_at_failed_coarse_and_fine_nodes(analysis):
    t,p=analysis['rows']
    assert (t['original_completed_SCFs_at_selected_coarse_nodes'],t['original_completed_SCFs_at_remaining_fine_nodes'])==(47,71)
    assert (p['original_completed_SCFs_at_selected_coarse_nodes'],p['original_completed_SCFs_at_remaining_fine_nodes'])==(58,114)
    # One coarse-corner hole stopped at10, not the12attempt ceiling.
    assert p['original_SCFs_at_missing_coarse_nodes']==46 and p['original_SCFs_at_fine_holes']==48
    for r in (t,p):
        assert r['original_completed_SCFs_at_selected_coarse_nodes']+r['original_completed_SCFs_at_remaining_fine_nodes']==r['source_candidate_SCFs']
    assert analysis['no_equal_budget_or_DFT_call_saving_or_B5_novelty_claim']


def test_same_cohort_signed_bias_and_empty_statistics_are_not_cherry_picked(analysis):
    t,p=analysis['rows']
    for m in a.MODELS:assert t['common_all_model_statistics'][m]['ids']==t['common_all_models_ids']
    assert t['common_all_model_statistics']['tensor_quadratic']['max_absolute_seen_origin_shape_error_meV_cell']==pytest.approx(2.003467715780971,abs=1e-11)
    assert t['model_origin_bias_meV_cell']['refitted_total_quadratic']>0
    assert t['common_all_model_statistics']['refitted_total_quadratic']['max_absolute_raw_error_meV_cell']!=t['common_all_model_statistics']['refitted_total_quadratic']['max_absolute_seen_origin_shape_error_meV_cell']
    for m in a.MODELS:
        stats=p['common_all_model_statistics'][m]
        assert stats['count']==0 and stats['max_absolute_seen_origin_shape_error_meV_cell'] is None
    assert analysis['all_labels_already_seen_before_this_exploratory_split']
    assert analysis['no_independent_prediction_or_branch_stability_or_physical_error_bound']
    assert analysis['new_DFT_scheduler_writes_reservations_or_independent_slots_used']==0
