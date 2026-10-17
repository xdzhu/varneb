"""Polynomial and full-gradient least-squares checks; no DFT or test labels."""
from pathlib import Path
import importlib.util
import json
import numpy as np
import pytest

PATH = Path(__file__).resolve().parents[1]/'benchmarks/hfo2_channels/20261008/G3_anharmonic_model_E265_20261018/fit.py'
spec = importlib.util.spec_from_file_location('E265_finite_fit_tests', PATH)
fit = importlib.util.module_from_spec(spec); spec.loader.exec_module(fit)


def test_polynomial_gradients_at_zero_and_general_points():
    exponents = np.array([[2,0,0],[1,1,0],[0,0,4],[1,0,2],[0,1,2]])
    points = np.array([[0.,0.,0.],[.2,-.3,.4]])
    _, jacobian = fit.monomials(points, exponents)
    for j in range(3):
        step = np.zeros(3); step[j] = 1e-6
        plus = fit.monomials(points+step, exponents)[0]
        minus = fit.monomials(points-step, exponents)[0]
        np.testing.assert_allclose(jacobian[:,j], (plus-minus)/2e-6, atol=2e-10, rtol=0)


def test_identified_energy_gradient_fit_recovers_known_polynomial():
    x = np.random.default_rng(41).normal(size=(30,3))
    exponents = np.array([[2,0,0],[1,1,0],[0,2,0],[0,0,3],[0,0,4]])
    values,jacobian = fit.monomials(x,exponents); known = np.array([.2,-.7,.3,.4,.02])
    energy = values@known; gradient = np.einsum('ndp,p->nd',jacobian,known)
    result = fit.solve(values,jacobian,energy,gradient,np.ones(len(x),dtype=bool))
    assert result['identified'] and result['rank']==5
    np.testing.assert_allclose(result['coefficients'],known,atol=2e-14,rtol=0)


def test_unidentified_fit_refuses_forecast_instead_of_minimum_norm_coefficients():
    values,jacobian = fit.monomials(np.zeros((4,3)),np.array([[2,0,0],[0,0,4]]))
    result = fit.solve(values,jacobian,np.zeros(4),np.zeros((4,3)),np.ones(4,dtype=bool))
    assert result['rank']==0 and not result['identified'] and result['coefficients'] is None


def test_exact_registered_sizes_no_degree_search_or_extra_bath_quadratic():
    harmonic,anharmonic=fit.powers(False),fit.powers(True)
    assert harmonic.shape==(3,36) and anharmonic.shape==(109,36)
    assert len(set(map(tuple,anharmonic)))==109
    assert np.all(anharmonic[:,3:].sum(axis=1)<=1)
    with pytest.raises(ValueError):fit.monomials(np.ones((2,3)),np.array([[-1,0,0]]))


@pytest.fixture(scope='module')
def trained():
    return json.loads((PATH.parent/'training_fit.json').read_bytes())


def test_actual_training_keeps_both_steps_and_worse_group_validation(trained):
    assert trained['native_training_observations']==81 and len(trained['models'])==4
    for step in (.01,.02):
        harmonic,anharmonic=[m for m in trained['models'] if m['step_A']==step]
        assert harmonic['fit']['rank']==3 and anharmonic['fit']['rank']==109
        assert anharmonic['in_sample']['energy_max_error_meV_fu']<harmonic['in_sample']['energy_max_error_meV_fu']
        assert anharmonic['leave_one_training_Q_family_out'][0]['metrics']['energy_max_error_meV_fu']>harmonic['leave_one_training_Q_family_out'][0]['metrics']['energy_max_error_meV_fu']
        assert [f['held_out_training_Q_id'] for f in anharmonic['leave_one_training_Q_family_out']]==list(range(5))
    assert not trained['independent_labels_seen'] and not trained['stable_DFT_sheet_TS_or_predictive_gain_certified']


def test_unclipped_negative_free_curvature_and_positive_bath_remain_distinct(trained):
    for model in trained['models']:
        transform=np.asarray(model['chart_transform_36_by_36']);extension=np.asarray(model['symmetric_measured_extension_eV_A2'])
        assert model['raw_R33_minimum_eV_A2']>0
        assert np.linalg.eigvalsh(transform[:,2:].T@extension@transform[:,2:])[0]<0
        np.testing.assert_allclose(transform[:,:2].T@extension@transform[:,:2],0,atol=2e-12,rtol=0)


def freezer_module():
    spec=importlib.util.spec_from_file_location('E265_frozen_evaluator_tests',PATH.parent/'freeze_P8.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def test_frozen_energy_dual_gradient_and_fixed_origin_match_finite_differences(trained):
    freezer=freezer_module();rng=np.random.default_rng(182)
    for model in trained['models']:
        origin=np.asarray(model['origin_internal_A']);e,g=freezer.evaluate(model,origin,fit)
        assert e==model['origin_energy_eV_cell']
        np.testing.assert_array_equal(g,model['origin_gradient_internal_eV_A'])
        x=origin+rng.normal(size=36)*.003;_,gradient=freezer.evaluate(model,x,fit)
        for j in (0,7,18,35):
            step=np.zeros(36);step[j]=2e-5
            derivative=(freezer.evaluate(model,x+step,fit)[0]-freezer.evaluate(model,x-step,fit)[0])/4e-5
            assert derivative==pytest.approx(gradient[j],abs=3e-7)


def test_complete_prospective_panel_preserves_refusal_and_no_independent_labels():
    freezer=freezer_module();result=freezer.freeze()
    saved=json.loads((PATH.parent/'P8_forecasts.json').read_bytes())
    assert result==saved and len(result['new_models'])==4
    assert [len(m['predictions']) for m in result['new_models']]==[8]*4
    all_points=[p for m in result['new_models'] for p in m['predictions']]
    assert sum(p['prediction_raw_eV_cell'] is not None for p in all_points)==31
    assert all(p['independent_DFT_label'] is None and p['model_unique_oracle_calls']<=12 for p in all_points)
    assert all(p['prediction_eV_fu']==p['prediction_raw_eV_cell']/4 for p in all_points if p['prediction_raw_eV_cell'] is not None)
    assert any(p['prediction_raw_eV_cell']<0 for p in all_points if p['prediction_raw_eV_cell'] is not None)
    assert len(result['original_reference_preserving_records'])==14
    assert len(result['original_standard_preserving_record']['targets'][0]['models'])==5
    assert result['original_conditional_2D_refusal_preserved']
    assert not result['stable_conditional_model_admission'] and not result['physical_B0_B5_barrier_panel_frozen']
    assert result['domain']['stable_conditional_predictive_domain']==[]
    assert result['domain']['maximum_total_internal_radius_A']==.5
    assert result['domain']['failed_training_candidates_preserved']==[0,1,6]
