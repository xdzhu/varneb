"""Analytic/injected guard records are not material SCFs."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest
from ase import Atoms

ROOT=Path(__file__).resolve().parents[1]
CASE=ROOT/'benchmarks/hfo2_channels/20261008/G3_bound_diagnostic_E168_20261015'
spec=importlib.util.spec_from_file_location('hfo2_bound_case',CASE/'diagnostics.py')
d=importlib.util.module_from_spec(spec);spec.loader.exec_module(d)
aspec=importlib.util.spec_from_file_location('hfo2_bound_arithmetic_case',CASE/'analyze.py')
a=importlib.util.module_from_spec(aspec);aspec.loader.exec_module(a)
SETTINGS=dict(maximum_internal_displacement_A=.5,minimum_distance_A=1.6,gradient_tolerance_eV_A=.03)


class Chart:
    def __init__(self,distance=2.):self.distance=distance;self.calls=0
    def to_atoms(self,x):
        self.calls+=1
        return Atoms('HH',positions=[[0,0,0],[self.distance,0,0]],cell=np.eye(3)*10,pbc=True)


@pytest.mark.parametrize('kind',('nonfinite','radius','accepted','distance','chart'))
def test_classifier_preserves_short_circuit_and_reason(kind):
    chart=Chart();x=np.zeros(36);basis=np.eye(40)[:,:36]
    if kind=='nonfinite':x[0]=np.nan
    if kind=='radius':x[0]=.50001
    if kind=='distance':chart.distance=1.
    if kind=='chart':
        def invalid(x):raise ValueError('explicit toy cell')
        chart.to_atoms=invalid
    got=d.guard_reason(chart,basis,x,SETTINGS)
    assert got['reason']==dict(nonfinite='nonfinite',radius='original_internal_radius',
        accepted='accepted',distance='minimum_distance',chart='chart_or_cell_invalid')[kind]
    assert got['safe']==(kind=='accepted')
    if kind in ('nonfinite','radius'):assert chart.calls==0 and got['minimum_distance_A'] is None


def test_wrong_dimension_is_not_classified_as_a_material_guard_failure():
    with pytest.raises(ValueError,match='36D'):d.guard_reason(Chart(),np.eye(40)[:,:36],np.zeros(35),SETTINGS)


@pytest.mark.parametrize('sign',[-1,1])
def test_trial_radial_diagnostic_is_not_an_optimizer_acceptance_or_KKT_certificate(sign):
    basis=np.eye(40)[:,:36];axes=np.eye(36)[:,:2];x=np.zeros(40);x[2]=.49
    g=np.zeros(40);g[2]=sign*.1
    record=dict(coordinates=x.tolist(),gradient=g.tolist(),enthalpy_eV=-5.)
    row=d.trial_rows(dict(basis=basis,axes=axes,center_energy=-4.),[record],.5,.03)[0]
    assert row['radius_margin_A']==pytest.approx(.01)
    assert row['radius_squared_steepest_descent_derivative_eV']==pytest.approx(-sign*.098)
    assert row['steepest_descent_points_outward']==(sign<0)
    assert not row['below_gradient_target'] and not row['accepted_iterate_or_final_candidate_inferred']
    assert row['energy_change_from_original_center_meV_cell']==-1000.


def test_missing_trials_are_not_replaced_by_a_final_candidate():
    with pytest.raises(ValueError,match='at least one'):
        d.trial_rows(dict(basis=np.eye(40)[:,:36],axes=np.eye(36)[:,:2]),[],.5,.03)


def replay_fixture():
    chart=Chart();basis=np.eye(40)[:,:36];axes=np.eye(36)[:,:2]
    x=np.zeros(40);x[2]=.49;g=np.zeros(40);g[2]=-.1
    records=[dict(coordinates=x.tolist(),gradient=g.tolist(),enthalpy_eV=-5.)]
    data=dict(chart=chart,basis=basis,axes=axes,center_energy=-4.)
    original=lambda chart,basis,x:d.guard_reason(chart,basis,x,SETTINGS)['safe']
    canary=SimpleNamespace(geometry_guard=original,SETTINGS=SETTINGS)
    audit=SimpleNamespace(worker=SimpleNamespace(canary=canary))
    row=dict(id=13,q_A=[0.,0.],status='optimizer_not_converged',reason='explicit toy budget',
             point=None,geometry_screens=dict(checks=2,rejected=1))
    def replay(data,row,records,start):
        assert canary.geometry_guard(chart,basis,np.zeros(36))
        bad=np.zeros(36);bad[2]=.6
        assert not canary.geometry_guard(chart,basis,bad)
        return dict(original_order_cache_only_replay_passed=True,new_DFT=0)
    audit.replay_node=replay
    return audit,data,row,records,original


def test_observer_preserves_and_restores_original_guard_and_does_not_promote_trial():
    audit,data,row,records,original=replay_fixture()
    got=d.observe_replay(audit,data,row,records,np.zeros(36))
    assert audit.worker.canary.geometry_guard is original
    assert got['guard_counts']==dict(accepted=1,original_internal_radius=1)
    assert got['reported_final_candidate'] is None and got['new_DFT_calls']==0
    assert not got['unconstrained_minimum_unreachable_or_boundary_KKT_or_branch_stability_certified']


def test_replay_failure_restores_the_live_imported_function():
    audit,data,row,records,original=replay_fixture()
    def bad(*args):raise RuntimeError('explicit toy mismatch')
    audit.replay_node=bad
    with pytest.raises(RuntimeError,match='mismatch'):d.observe_replay(audit,data,row,records,np.zeros(36))
    assert audit.worker.canary.geometry_guard is original


def test_classifier_mismatch_is_not_silently_used_as_a_repaired_guard():
    audit,data,row,records,_=replay_fixture()
    original=lambda *args:False
    audit.worker.canary.geometry_guard=original
    with pytest.raises(ValueError,match='classifier differs'):d.observe_replay(audit,data,row,records,np.zeros(36))
    assert audit.worker.canary.geometry_guard is original


def test_actual_native_prefix_keeps_all_twelve_trials_and_exact_radius_rejections():
    result=a.analyze()
    assert result['original_prefix_nodes']==14 and result['original_prefix_candidates']==13
    assert result['raw_SCFs_reparsed_by_native_driver']==60
    assert len(result['holes'])==1
    hole=result['holes'][0]
    assert hole['point_id']==13 and hole['q_A']==[-.2,-.1]
    assert hole['guard_counts']==dict(accepted=13,original_internal_radius=35)
    assert len(hole['trial_rows'])==12 and hole['no_completed_trial_below_target']
    assert hole['minimum_observed_trial_gradient_norm_eV_A']>.03
    assert not hole['boundary_KKT_or_unreachable_unconstrained_minimum_or_unique_cause_certified']
    assert not result['original_input_domain_force_target_and_runtime_changed']
    assert not result['physical_reserved_geometry_generated_or_read']


def test_actual_last_trial_still_has_a_large_constant_radius_tangent_gradient():
    row=a.analyze()['holes'][0]['last_completed_trial']
    assert row['internal_norm_A']==pytest.approx(.4998853734687307,abs=1e-14)
    assert row['released_gradient_norm_eV_A']==pytest.approx(.07359637132369291,abs=1e-14)
    assert row['radial_gradient_component_eV_A']==pytest.approx(-.013612340116481767,abs=1e-14)
    assert row['constant_radius_tangent_gradient_norm_eV_A']==pytest.approx(.07232655161535154,abs=1e-14)
    assert row['steepest_descent_points_outward'] and not row['accepted_iterate_or_final_candidate_inferred']


def test_scalar_split_rejects_a_false_physical_gradient_identity():
    with pytest.raises(ValueError,match='exceeds total'):a.radial_split(1.,.1,1.)


@pytest.mark.parametrize('values',[(0.,1.,0.),(-1.,1.,0.),(1.,-1.,0.),(1.,1.,np.nan)])
def test_scalar_split_rejects_undefined_or_nonfinite_data(values):
    with pytest.raises(ValueError,match='finite positive'):a.radial_split(*values)


def test_arithmetic_CLI_never_overwrites_or_reparses_before_refusing(tmp_path,monkeypatch):
    path=tmp_path/'old.json';path.write_text('preserve')
    monkeypatch.setattr('sys.argv',['diagnose','--output',str(path)])
    monkeypatch.setattr(a,'analyze',lambda:pytest.fail('must refuse before reading source'))
    with pytest.raises(FileExistsError):a.main()
    assert path.read_text()=='preserve'


def test_case_table_preserves_the_first_minimum_norm_and_last_trials_without_candidate_promotion():
    hole=a.analyze()['holes'][0];rows=hole['trial_rows']
    chosen=[rows[0],min(rows,key=lambda r:r['released_gradient_norm_eV_A']),rows[-1]]
    assert [r['trial_index'] for r in chosen]==[0,5,11]
    text=(CASE/'README.md').read_text(encoding='utf-8')
    for r in chosen:
        line=(f"| {r['trial_index']} | {r['internal_norm_A']:.6f} | "
              f"{r['released_gradient_norm_eV_A']:.6f} | "
              f"{r['constant_radius_tangent_gradient_norm_eV_A']:.6f} | "
              f"{r['energy_change_from_original_center_meV_cell']:.6f} |")
        assert text.splitlines().count(line)==1
    assert 'No trial is promoted' in text and 'all12 remain' in text
