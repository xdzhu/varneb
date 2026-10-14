"""Analytic/mocked reader checks, not HfO2 material or TS certification."""
import copy
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from ase.calculators.singlepoint import SinglePointCalculator

from scripts import audit_hfo2_G3_curvature as audit


@pytest.fixture(scope='module')
def manifest():
    path=Path(__file__).resolve().parents[1]/'benchmarks/hfo2_channels/20261008/G3_curvature_E150_20261014/manifest.json'
    return json.loads(path.read_text())


def quadratic(manifest, n=148):
    basis=np.asarray(manifest['centers']['T']['directions'])
    rng=np.random.default_rng(28)
    rotation=np.linalg.qr(rng.normal(size=(36,36)))[0]
    internal=rotation@np.diag(np.r_[-2.,np.linspace(.4,5,35)])@rotation.T
    reduced=np.zeros((37,37)); reduced[:36,:36]=internal
    reduced[-1,-1]=-20.  # controlled coordinate must not add an internal instability
    reduced[-1,:36]=reduced[:36,-1]=rng.normal(size=36)*.03
    matrix=basis@reduced@basis.T
    gradient=rng.normal(size=40)*.01
    _,probes=audit.study.validate(manifest,'T')
    rows=[dict(status='complete',coordinates=p.delta_A.tolist(),
        enthalpy_eV=float(10.+gradient@p.delta_A+p.delta_A@matrix@p.delta_A/2),
        gradient=(gradient+matrix@p.delta_A).tolist()) for p in probes[:n]]
    return gradient, rows


@pytest.mark.parametrize('n',[0,1,2,29,73,74,75,76,147,148])
def test_incomplete_blocks_do_not_invent_missing_columns_or_full_index(manifest,n):
    gradient,rows=quadratic(manifest,n)
    result=audit.diagnostics(manifest,'T',10.,gradient,rows)
    assert result['completed_records']==n and result['planned_records']==148
    assert (result['two_step_comparison'] is not None)==(n==148)
    assert len(result['two_step_lines'])==max(0,min(74,n-74))//2
    for b in result['blocks']:
        k=b['complete_pairs']
        assert np.shape(b['full_hessian_action_eV_A2'])==(40,k)
        assert np.shape(b['raw_projected_eV_A2'])==(k,k)
        assert (b['fixed_strain_internal_eigenvalues_eV_A2'] is not None)==(k==37)
    assert not result['physical_error_bound_or_TS_certificate']
    assert not result['selected_conditional_plane']


def test_full_coupled_quadratic_matches_existing_energy_gradient_and_work_checks(manifest):
    gradient,rows=quadratic(manifest)
    result=audit.diagnostics(manifest,'T',10.,gradient,rows)
    comparison=result['two_step_comparison']
    assert comparison['internal_counts_below_zero']==[1,1]
    assert comparison['lowest_internal_direction_absolute_overlap']==pytest.approx(1.)
    assert comparison['internal_operator_step_change_eV_A2']<1e-12
    for b in result['blocks']:
        assert np.max(np.abs(b['energy_gradient_diagonal_defect_eV_A2']))<1e-10
    for line in result['two_step_lines']:
        assert np.max(np.abs(line['five_point_work']['energy_minus_gradient_work']))<1e-13
    assert result['two_step_lines'][-1]['control_direction']
    json.dumps(result,allow_nan=False)


@pytest.mark.parametrize('fault',['reorder','coordinate','gradient','energy','extra'])
def test_malformed_complete_prefix_is_not_repaired(manifest,fault):
    gradient,rows=quadratic(manifest,4)
    if fault=='reorder': rows[0],rows[1]=rows[1],rows[0]
    if fault=='coordinate': rows[0]['coordinates'][0]+=.01
    if fault=='gradient': rows[0]['gradient'][0]=float('nan')
    if fault=='energy': rows[0]['enthalpy_eV']=float('inf')
    if fault=='extra': rows*=38
    with pytest.raises(ValueError): audit.diagnostics(manifest,'T',10.,gradient,rows)


def reader_fixture(tmp_path,manifest,monkeypatch,n=3):
    chart,probes=audit.study.validate(manifest,'T')
    monkeypatch.setattr(audit.study,'native_preflight',lambda *args:(chart,probes,0.,np.zeros(40)))
    work=tmp_path/'live'; directory=work/'evaluations'; directory.mkdir(parents=True)
    (work/'center.json').write_text(json.dumps(dict(cached_energy_eV_cell=0.,
        cached_full_gradient_eV_A=np.zeros(40).tolist(),new_center_SCFs=0,job_id='MOCK-NOT-SLURM')))
    contract=dict(coordinate_convention=chart.coordinate_convention,mechanical_boundary=chart.mechanical_contract,
        pressure_eV_per_A3=0.,calculator_id='original-six-HfO2-PBE-100Ry-full10auDZP-Gamma2x2x2',
        validation_id='E150-native-six-input-DSIZE32-EFS-one-SCF-v1')
    digest=hashlib.sha256(json.dumps(contract,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    contract['contract_sha256']=digest
    (directory/'contract.json').write_text(json.dumps(contract))
    validator=dict(explicit_analytic_mock_not_native=True,elapsed_seconds=1.)
    raw=dict(energy=0.,forces=np.zeros((12,3)),stress=np.zeros(6))
    monkeypatch.setattr(audit.study,'audited_results',lambda p:copy.deepcopy(raw))
    validated=[]
    monkeypatch.setattr(audit.study,'result_validator',lambda i,a,p:(validated.append(i) or validator.copy()))
    paths=[]
    for i,p in enumerate(probes[:n]):
        d=hashlib.sha256(p.delta_A.tobytes()).hexdigest()
        point=directory/f'eval-{i:06d}-{d[:12]}'; (point/'scf_000000').mkdir(parents=True)
        a=p.atoms; a.calc=SinglePointCalculator(a,**raw)
        e,g=chart.energy_and_gradient(p.delta_A,a)
        path=point/'result.json'; paths.append(path)
        path.write_text(json.dumps(dict(status='complete',coordinate_sha256=d,contract_sha256=digest,
            coordinates=p.delta_A.tolist(),enthalpy_eV=e,gradient=g.tolist(),
            stress_eV_per_A3=a.get_stress(voigt=False).tolist(),validation=validator)))
    return work,paths,validated


def test_actual_reader_flow_reparses_each_completed_point_with_explicit_raw_mocks(tmp_path,manifest,monkeypatch):
    work,paths,validated=reader_fixture(tmp_path,manifest,monkeypatch)
    result=audit.audit(manifest,'T',work)
    assert validated==[0,1,2] and result['raw_SCFs_reparsed']==3
    assert result['observed_native_SCF_directories']==3
    assert result['sum_audited_call_seconds']==3.
    assert result['new_DFT_or_scheduler_operations']==0
    assert result['allocation_core_hours'] is None
    assert not result['physical_reserved_geometry_generated_or_read']


@pytest.mark.parametrize('fault',['hole','cached_gradient','cached_energy','cached_stress','cached_validation','contract','failure'])
def test_reader_refuses_or_exposes_faults_without_touching_native_files(tmp_path,manifest,monkeypatch,fault):
    work,paths,_=reader_fixture(tmp_path,manifest,monkeypatch)
    if fault=='hole': paths[0].unlink()
    elif fault=='contract':
        path=work/'evaluations/contract.json'; c=json.loads(path.read_text());c['pressure_eV_per_A3']=1.
        path.write_text(json.dumps(c))
    elif fault=='failure': (paths[-1].parent/'failure.json').write_text('{"explicit_mock":true}')
    else:
        row=json.loads(paths[0].read_text())
        if fault=='cached_gradient': row['gradient'][0]+=.1
        if fault=='cached_energy': row['enthalpy_eV']+=.1
        if fault=='cached_stress': row['stress_eV_per_A3'][0][0]+=.1
        if fault=='cached_validation': row['validation']['elapsed_seconds']=2.
        paths[0].write_text(json.dumps(row))
    before={p:str(p.read_bytes()) for p in work.rglob('*') if p.is_file()}
    if fault=='failure':
        assert audit.audit(manifest,'T',work)['status']=='native_evaluation_failure_present'
    else:
        with pytest.raises(ValueError): audit.audit(manifest,'T',work)
    assert {p:str(p.read_bytes()) for p in work.rglob('*') if p.is_file()}==before
