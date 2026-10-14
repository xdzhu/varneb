"""Block-curvature arithmetic is not a material stability or TS certificate."""
from copy import deepcopy
import hashlib
import json

import numpy as np
import pytest

from scripts import analyze_hfo2_G3_atomic_cell as analysis


@pytest.fixture
def example():
    # Explicit toy embedding, not a physical translation gauge or DFT record.
    b=np.eye(40)[:,[*range(33),36,37,38]]
    h=np.diag(np.r_[np.ones(33),[2.,3.,4.]])
    g=np.zeros(40)
    return h,b,g


def test_zero_coupling_does_not_soften_atomic_curvature(example):
    h,b,g=example;h[0,0]=-.3
    result=analysis.raw_atomic_cell_signature(h,b,g)
    formal=result['formal_cell_stationary_quadratic']
    assert result['atomic_raw_negative_count']==result['internal_raw_negative_count']==1
    np.testing.assert_allclose(result['atomic_fixed_cell_eigenvalues_eV_A2'],
        formal['atomic_curvature_eigenvalues_eV_A2'],atol=2e-12,rtol=0)
    assert result['negative_subspace']['open_cell_projector_trace']==0.
    assert result['negative_subspace']['cross_Rayleigh_trace_eV_A2']==0.


def test_coupling_can_create_an_additional_negative_direction_not_hardcoded(example):
    h,b,g=example;h[0,0]=-1.;h[33,33]=1.;h[1,33]=h[33,1]=2.
    result=analysis.raw_atomic_cell_signature(h,b,g)
    assert result['atomic_raw_negative_count']==1
    assert result['internal_raw_negative_count']==2
    assert result['formal_cell_stationary_quadratic']['atomic_raw_negative_count']==2
    assert result['formal_cell_stationary_quadratic']['atomic_curvature_eigenvalues_eV_A2'][0]==pytest.approx(-3.)
    trace=result['negative_subspace']
    assert trace['open_cell_projector_trace']==pytest.approx(.5)
    assert (trace['atomic_Rayleigh_trace_eV_A2']+trace['open_cell_Rayleigh_trace_eV_A2']+
            trace['cross_Rayleigh_trace_eV_A2'])==pytest.approx(trace['sum_negative_eigenvalues_eV_A2'])


def test_nonstationary_linear_cell_term_is_retained(example):
    h,b,g=example;h[0,33]=h[33,0]=.5;g[0]=.1;g[36]=.3
    formal=analysis.raw_atomic_cell_signature(h,b,g)['formal_cell_stationary_quadratic']
    assert formal['nonzero_gradient_cell_offset_norm_A']==pytest.approx(.15)
    assert formal['constant_energy_offset_eV_cell']==pytest.approx(-.0225)
    assert formal['effective_atomic_gradient_norm_eV_A']==pytest.approx(.025)
    assert not formal['eliminated_gradient_assumed_zero']
    assert formal['raw_quadratic_only_not_DFT_prediction_or_measured_error_bound']


@pytest.mark.parametrize('cell_curvature',[-1.,0.])
def test_nonpositive_cell_block_is_never_clipped_or_pseudoinverted(example,cell_curvature):
    h,b,g=example;h[33,33]=cell_curvature
    result=analysis.raw_atomic_cell_signature(h,b,g)
    assert result['formal_cell_stationary_quadratic'] is None
    assert result['open_cell_raw_nonpositive_count']==1


def compare_invariants(first,second):
    for key in ('internal_eigenvalues_eV_A2','atomic_fixed_cell_eigenvalues_eV_A2','open_cell_eigenvalues_eV_A2'):
        np.testing.assert_allclose(first[key],second[key],atol=2e-12,rtol=0)
    for key in ('internal_raw_negative_count','atomic_raw_negative_count','open_cell_raw_nonpositive_count'):
        assert first[key]==second[key]
    assert second['atomic_cell_coupling_operator_norm_eV_A2']==pytest.approx(
        first['atomic_cell_coupling_operator_norm_eV_A2'],abs=2e-12,rel=0)
    for key,value in first['negative_subspace'].items():
        assert second['negative_subspace'][key]==pytest.approx(value,abs=2e-12,rel=0)
    if first['formal_cell_stationary_quadratic'] is None:
        assert second['formal_cell_stationary_quadratic'] is None
    else:
        for key,value in first['formal_cell_stationary_quadratic'].items():
            if isinstance(value,bool):assert second['formal_cell_stationary_quadratic'][key] is value
            else:np.testing.assert_allclose(second['formal_cell_stationary_quadratic'][key],value,atol=2e-12,rtol=0)


@pytest.mark.parametrize('degenerate',[False,True])
def test_internal_basis_rotation_does_not_change_physical_partition_traces(example,degenerate):
    h,b,g=example
    h[0,0]=-1.
    if degenerate:h[33,33]=-1.
    else:h[1,33]=h[33,1]=2.;g[36]=.2
    q,_=np.linalg.qr(np.random.default_rng(381).normal(size=(36,36)))
    first=analysis.raw_atomic_cell_signature(h,b,g)
    second=analysis.raw_atomic_cell_signature(q.T@h@q,b@q,g)
    compare_invariants(first,second)


@pytest.mark.parametrize('fault',['asymmetric','nonfinite','complex','wrong_h_shape','wrong_g_shape',
                                'nonorthogonal','control_in_basis','wrong_partition'])
def test_incomplete_or_wrong_coordinate_problem_is_refused(example,fault):
    h,b,g=deepcopy(example)
    if fault=='asymmetric':h[0,1]=.2
    elif fault=='nonfinite':h[0,0]=float('nan')
    elif fault=='complex':h=h.astype(complex);h[0,0]+=1j
    elif fault=='wrong_h_shape':h=h[:35,:35]
    elif fault=='wrong_g_shape':g=g[:39]
    elif fault=='nonorthogonal':b[:,0]*=2
    elif fault=='control_in_basis':b[36,33]=0.;b[39,33]=1.
    else:b=np.eye(40)[:,:36]
    with pytest.raises(ValueError):analysis.raw_atomic_cell_signature(h,b,g)


@pytest.mark.parametrize('channel,step,atomic,joint,cell_trace',[
    ('T',.01,-.24016859376397562,-.27619342870577585,.020577327206936582),
    ('T',.02,-.23780578381134732,-.2735964758435756,.02061046493973889),
    ('preserving',.01,-1.271818022831291,-1.2718188494432225,1.0402204001767562e-6),
    ('preserving',.02,-1.2818154704685132,-1.2818205410155206,3.092475648182485e-6)])
def test_complete_actual_two_step_blocks_reproduce_signature_and_boundaries(channel,step,atomic,joint,cell_trace):
    result=analysis.analyze()
    row=next(r for r in result['rows'] if (r['channel'],r['step_A'])==(channel,step))
    assert row['atomic_fixed_cell_eigenvalues_eV_A2'][0]==pytest.approx(atomic,abs=2e-12,rel=0)
    assert row['internal_eigenvalues_eV_A2'][0]==pytest.approx(joint,abs=2e-12,rel=0)
    assert row['negative_subspace']['open_cell_projector_trace']==pytest.approx(cell_trace,abs=2e-12,rel=0)
    assert row['atomic_raw_negative_count']==row['internal_raw_negative_count']==2
    assert row['open_cell_raw_nonpositive_count']==0
    assert row['formal_cell_stationary_quadratic']['atomic_raw_negative_count']==2
    assert row['algebraic_zero_signs_not_physical_resolution']
    assert not row['phonon_identity_stationarity_branch_stability_or_barrier_decomposition_certified']
    assert result['new_DFT_or_scheduler_operations']==0 and not result['new_raw_SCF_reparse']
    assert not result['physical_reserved_geometry_generated_or_read']
    assert not result['conditional_axes_final_representatives_or_forecasts_selected']
    assert result['cell_scale_A']==5.12968067458423


def test_existing_output_refuses_before_any_data_read(tmp_path,monkeypatch):
    path=tmp_path/'existing.json';path.write_text('preserve')
    monkeypatch.setattr('sys.argv',['analysis','--output',str(path)])
    monkeypatch.setattr(analysis,'analyze',lambda:pytest.fail('must refuse before reading data'))
    with pytest.raises(FileExistsError):analysis.main()
    assert path.read_text()=='preserve'


def test_saved_analysis_contains_no_individual_eigenvectors_or_author_arrays():
    text=json.dumps(analysis.analyze())
    assert 'eigenvectors' not in text and 'source_author' not in text


def test_frozen_actual_report_replays_every_field_without_DFT():
    case=analysis.BASE/'G3_atomic_cell_signature_E164_20261015'
    raw=(case/'analysis.json').read_bytes()
    assert hashlib.sha256(raw).hexdigest()=='ef7798df3e71fbbc3e86bd468b608e4bf732cb13cd048e1518b8608d4753db3a'
    def same(actual,expected):
        if isinstance(expected,dict):
            assert actual.keys()==expected.keys()
            for key in expected:same(actual[key],expected[key])
        elif isinstance(expected,list):
            assert len(actual)==len(expected)
            for a,e in zip(actual,expected):same(a,e)
        elif isinstance(expected,float):assert actual==pytest.approx(expected,abs=2e-12,rel=0)
        else:assert type(actual) is type(expected) and actual==expected
    same(analysis.analyze(),json.loads(raw))


def test_case_table_is_rounded_from_all_four_actual_center_step_rows():
    case=analysis.BASE/'G3_atomic_cell_signature_E164_20261015'
    report=json.loads((case/'analysis.json').read_bytes())
    text=(case/'README.md').read_text(encoding='utf-8')
    for row in report['rows']:
        label={'T':'T','preserving':'Preserving'}[row['channel']]
        expected=(f"| {label} | {row['step_A']:.2f} | {row['atomic_fixed_cell_eigenvalues_eV_A2'][0]:.6f} | "
            f"{row['internal_eigenvalues_eV_A2'][0]:.6f} | {row['open_cell_eigenvalues_eV_A2'][0]:.6f} | "
            f"{row['negative_subspace']['open_cell_projector_trace']:.6g} |")
        assert text.splitlines().count(expected)==1
