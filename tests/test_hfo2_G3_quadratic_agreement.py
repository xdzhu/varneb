"""Retrospective agreement is not independent prediction or an error bound."""
import importlib.util
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from vcneb.quadratic_reduction import condition_quadratic_energy

ROOT=Path(__file__).resolve().parents[1]
CASE=ROOT/'benchmarks/hfo2_channels/20261008/G3_quadratic_agreement_E165_20261015'
spec=importlib.util.spec_from_file_location('hfo2_G3_quadratic_agreement_case',CASE/'analyze.py')
analysis=importlib.util.module_from_spec(spec);spec.loader.exec_module(analysis)


@pytest.fixture(scope='module')
def actual():return analysis.analyze()


@pytest.fixture
def exact_model_and_shifted_observations():
    h=np.diag([-1.,-2.,1.]);q=np.eye(3)[:,:2];r=np.eye(3)[:,2:]
    model=condition_quadratic_energy(h,np.array([.2,-.1,.3]),-4.,q,r,stability_floor=0.)
    qs=[[0.,0.],[-.1,0.],[0.,-.1],[0.,.1],[.1,0.]]
    points=[dict(id=i,q_A=x,energy_eV_cell=-4.+model.energy_change(x)+.01,
        retained_gradient_eV_A=(model.retained_gradient(x)+[.02,0.]).tolist()) for i,x in enumerate(qs)]
    return model,points


def test_seen_origin_alignment_hides_constant_error_but_not_gradient_error(exact_model_and_shifted_observations):
    model,points=exact_model_and_shifted_observations
    rows=analysis.evaluate_points(model,points,-4.)
    np.testing.assert_allclose([r['energy_error_meV_cell'] for r in rows],-10.,atol=2e-12,rtol=0)
    np.testing.assert_allclose([r['zero_aligned_shape_error_meV_cell'] for r in rows],0.,atol=2e-12,rtol=0)
    np.testing.assert_allclose([r['retained_partial_gradient_error_norm_eV_A'] for r in rows],.02,atol=2e-12,rtol=0)


def test_axial_energy_and_gradient_secants_reproduce_only_the_toy_quadratic(exact_model_and_shifted_observations):
    _,points=exact_model_and_shifted_observations
    result=analysis.axial_differences(points)
    np.testing.assert_allclose(result['candidate_energy_axial_second_differences_eV_A2'],[-1.,-2.],atol=2e-12,rtol=0)
    np.testing.assert_allclose(result['raw_retained_partial_gradient_secants_eV_A2'],np.diag([-1.,-2.]),atol=2e-12,rtol=0)
    assert not result['infinitesimal_conditional_Hessian_or_physical_error_bound']


def test_common_statistics_refuse_missing_declared_points(exact_model_and_shifted_observations):
    model,points=exact_model_and_shifted_observations
    rows=analysis.evaluate_points(model,points,-4.)
    with pytest.raises(ValueError,match='common numeric-Q'):analysis.statistics(rows,{0,1,5})


@pytest.mark.parametrize('channel,step,offset,max_shape',[
    ('T',.01,-.8669582250244786,8.524106080376756),
    ('T',.02,-.8744289915789629,8.58780036999137),
    ('preserving',.01,-54.52122170979342,20.52183879154757),
    ('preserving',.02,-54.70582849920575,20.72063953107716)])
def test_actual_common_cohort_keeps_nonzero_origin_bias_and_shape_disagreement(actual,channel,step,offset,max_shape):
    row=next(r for r in actual['rows'] if (r['channel'],r['curvature_probe_step_A'])==(channel,step))
    assert row['raw_formal_quadratic_origin_release_meV_cell']==pytest.approx(offset,abs=2e-10,rel=0)
    assert row['common_numeric_Q_statistics']['max_absolute_zero_aligned_shape_error_meV_cell']==pytest.approx(max_shape,abs=2e-10,rel=0)
    assert row['common_numeric_Q_statistics']['ids']==list(range(6))
    assert len(row['points'])=={'T':9,'preserving':6}[channel]
    assert row['points'][0]['energy_error_meV_cell']!=0.
    assert row['points'][0]['zero_aligned_shape_error_meV_cell']==0.
    if step==.01:assert row['formal_Schur_softening_Frobenius_norm_eV_A2']<1e-24


@pytest.mark.parametrize('channel,expected',[
    ('T',[-.40512781815778,-.24308684332936537]),
    ('preserving',[-.33776329619286116,-.18457582882547285])])
def test_actual_axial_candidate_differences_are_not_labelled_Hessians(actual,channel,expected):
    row=actual['candidate_axial_differences'][channel]
    np.testing.assert_allclose(row['candidate_energy_axial_second_differences_eV_A2'],expected,atol=2e-10,rtol=0)
    assert row['step_A']==.1 and not row['infinitesimal_conditional_Hessian_or_physical_error_bound']


def test_actual_source_cohort_and_boundaries_are_not_independent_success(actual):
    assert actual['available_nodes']=={'T':9,'preserving':6}
    assert actual['common_numeric_Q_ids']==list(range(6))
    assert actual['common_numeric_Q_does_not_mean_same_physical_configurations']
    assert actual['all_model_parameters_from_original_center_no_coefficient_fit']
    assert actual['zero_alignment_uses_seen_origin_and_is_only_a_shape_diagnostic']
    assert actual['algebraic_stability_floor_zero_not_physical_resolution']
    assert actual['new_DFT_or_scheduler_operations']==0 and not actual['new_raw_SCF_reparse']
    assert not actual['conditional_stability_full_faces_TS_or_material_prediction_gain_certified']
    assert not actual['physical_reserved_geometry_generated_or_read']


def test_case_table_uses_all_four_original_models_without_error_based_selection(actual):
    text=(CASE/'README.md').read_text(encoding='utf-8')
    for row in actual['rows']:
        label={'T':'T','preserving':'Preserving'}[row['channel']]
        expected=(f"| {label} | {row['curvature_probe_step_A']:.2f} | "
            f"{row['observed_DFT_candidate_origin_release_meV_cell']:.6f} | "
            f"{row['raw_formal_quadratic_origin_release_meV_cell']:.6f} | "
            f"{row['common_numeric_Q_statistics']['max_absolute_zero_aligned_shape_error_meV_cell']:.6f} |")
        assert text.splitlines().count(expected)==1


def test_existing_output_refuses_before_source_reads(tmp_path,monkeypatch):
    path=tmp_path/'old.json';path.write_text('preserve')
    monkeypatch.setattr('sys.argv',['agreement','--output',str(path)])
    monkeypatch.setattr(analysis,'analyze',lambda:pytest.fail('must refuse first'))
    with pytest.raises(FileExistsError):analysis.main()
    assert path.read_text()=='preserve'


def test_changed_native_prefix_pin_refuses_without_reanalysis(monkeypatch):
    monkeypatch.setattr(analysis,'PINS',dict(analysis.PINS,T='0'*64))
    with pytest.raises(ValueError,match='prefix bytes differ'):analysis.analyze()


def test_saved_complete_report_replays_all_fields_with_strict_source_bytes(actual):
    raw=(CASE/'analysis.json').read_bytes()
    assert hashlib.sha256(raw).hexdigest()=='c2e51c37b185a00676a6ff21660ea132384c0a22cfac2cebd18e30a5e5a9b745'
    def same(a,e):
        if isinstance(e,dict):
            assert a.keys()==e.keys()
            for key in e:same(a[key],e[key])
        elif isinstance(e,list):
            assert len(a)==len(e)
            for av,ev in zip(a,e):same(av,ev)
        elif isinstance(e,float):assert a==pytest.approx(e,abs=2e-10,rel=0)
        else:assert type(a) is type(e) and a==e
    same(json.loads(json.dumps(actual)),json.loads(raw))
