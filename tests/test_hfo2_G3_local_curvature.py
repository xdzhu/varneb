"""Finite native training probe admission; mocks never constitute HfO2 DFT."""
import copy
import json
from pathlib import Path

import numpy as np
import pytest
from ase.calculators.singlepoint import SinglePointCalculator

import scripts.hfo2_G3_local_curvature as study


@pytest.fixture(scope='module')
def manifest(tmp_path_factory):
    target=tmp_path_factory.mktemp('registered')/'manifest.json'
    return study.prepare(target)


def test_actual_complete_training_selection_is_not_a_physical_ordering_claim(manifest):
    assert set(manifest['centers']) == {'T','preserving'}
    assert manifest['ordinary_G2_paths'] == 8
    assert manifest['barrier_error_bounds_measured'] is False
    assert manifest['final_predictive_edge_selection_locked'] is False
    assert manifest['new_production_chains'] == 0
    assert manifest['physical_holdout_generated_or_read'] is False
    assert manifest['total_planned_raw_SCF_attempts'] == 296
    assert manifest['grid_nodes_used'] == manifest['independent_face_test_points_used'] == 0


@pytest.mark.parametrize('channel', ['T','preserving'])
def test_full_native_stencil_has_exact_control_internal_space_and_safe_geometry(manifest,channel):
    chart,probes=study.validate(manifest,channel)
    assert len(probes)==148 and chart.size==40
    assert chart.internal_relaxation_basis().shape==(40,36)
    for p in probes:
        strain=p.delta_A[-1]/chart.strain_scale_A
        assert not np.isclose(strain,.005,atol=1e-12,rtol=0)
        if p.direction_index==36:
            assert strain==pytest.approx(p.sign*p.step_A/chart.strain_scale_A)
        else:
            assert strain==0.
        np.testing.assert_allclose(p.atoms.cell[:2],chart.reference_cell[:2]*(1+strain),atol=1e-12)
        assert p.atoms.calc is None
        np.testing.assert_allclose(chart.from_atoms(p.atoms),p.delta_A,atol=1e-12)


@pytest.mark.parametrize('key,value', [('pressure_GPa',1.),('anchor_strain',.005),
    ('climb',True),('MPI_ranks_per_allocation',16),('ordinary_G2_paths',7),
    ('maximum_concurrent_allocations',3),('planned_raw_SCF_attempts_per_center',149),
    ('step_amplitudes_A',[.01,.03]),('total_planned_raw_SCF_attempts',True)])
def test_parameter_budget_or_holdout_drift_is_refused(manifest,key,value):
    altered=copy.deepcopy(manifest); altered[key]=value
    with pytest.raises(ValueError,match='contract differs'): study.validate(altered,'T')


@pytest.mark.parametrize('fault',['source','code','stencil','basis','physical'])
def test_source_geometry_basis_and_original_contract_cannot_be_replaced(manifest,fault):
    altered=copy.deepcopy(manifest)
    if fault=='source': altered['centers']['T']['source_job_id']='27719610'
    if fault=='code': altered['normalized_source_text_sha256'][study.CODE[0]]='0'*64
    if fault=='stencil': altered['centers']['T']['probe_coordinates_A'][0][0]+=.001
    if fault=='basis': altered['centers']['T']['directions'][-1][0]=.1
    if fault=='physical': altered['physical_inputs_sha256']['INPUT']='0'*64
    with pytest.raises(ValueError): study.validate(altered,'T')


def test_preparation_does_not_overwrite_or_accept_partial_training(tmp_path,monkeypatch):
    target=tmp_path/'existing.json'; target.write_text('original')
    with pytest.raises(FileExistsError): study.prepare(target)
    assert target.read_text()=='original'
    import scripts.audit_hfo2_clamped_network as network
    monkeypatch.setattr(network,'audit',lambda: {'response': {'sampled_minima': {'0.0': {
        'sampled_switching_candidate':'reversing','sampled_escape_candidate':'M'}}}})
    with pytest.raises(ValueError,match='training selection differs'): study.prepare(tmp_path/'different.json')


def test_login_or_wrong_resources_never_invoke_preflight_or_calculator(tmp_path,manifest,monkeypatch):
    calls=[]
    monkeypatch.setattr(study,'native_preflight',lambda *args: calls.append(True))
    monkeypatch.delenv('RUN_DFT',raising=False)
    with pytest.raises(ValueError,match='allocation'): study.run(manifest,'T',tmp_path/'never')
    assert calls==[] and not (tmp_path/'never').exists()


def test_static_driver_all_points_and_cost_with_explicit_analytic_mocks_not_DFT(tmp_path,manifest,monkeypatch):
    chart,probes=study.validate(manifest,'T')
    monkeypatch.setattr(study,'RUNS',tmp_path)
    monkeypatch.setattr(study,'native_preflight',lambda *args: (chart,probes,0.,np.zeros(40)))
    calls=[]
    def factory(**kwargs):
        assert kwargs['command'].startswith('mpirun -np 32 ')
        def calculate(index,atoms,directory):
            calls.append(index)
            return SinglePointCalculator(atoms,energy=0.,forces=np.zeros((12,3)),stress=np.zeros(6))
        return calculate
    monkeypatch.setattr(study,'make_factory',factory)
    monkeypatch.setattr(study,'result_validator',lambda *args: dict(elapsed_seconds=1.,explicit_analytic_mock=True))
    for name,value in dict(RUN_DFT='1',SLURM_JOB_PARTITION='hfacnormal01',SLURM_JOB_NUM_NODES='1',
        SLURM_NTASKS='1',SLURM_CPUS_PER_TASK='32',SLURM_JOB_ID='TEST-NOT-SLURM',
        OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1').items():
        monkeypatch.setenv(name,value)
    result=study.run(manifest,'T',tmp_path/'analytic-only')
    assert calls==list(range(148))
    assert result['sum_native_call_seconds']==148.  # injected analytic time, not measured DFT
    assert result['call_time_proxy_core_hours']==pytest.approx(148*32/3600)
    assert result['allocation_core_hours'] is None
    assert result['conditional_surface_or_predictions_completed'] is False
    for block in result['blocks']:
        np.testing.assert_array_equal(block['raw_projected_eV_A2'],np.zeros((37,37)))
        assert len(block['fixed_strain_internal_eigenvalues_eV_A2'])==36
    with pytest.raises(ValueError,match='fresh'): study.run(manifest,'T',tmp_path/'analytic-only')


@pytest.mark.parametrize('fault', ['none','extra_SCF','energy','force','stress','input','zero_time','nan_time','bool_time'])
def test_native_validator_checks_raw_EFS_physical_identity_and_call_count(tmp_path,monkeypatch,fault):
    import shutil
    from scripts.analyze_hfo2_clamped_reference import terminal_images
    base=study.ROOT/'benchmarks/hfo2_channels/20261008'
    case,pin,_=study.SOURCES['T']
    images,obs,_,receipt=terminal_images(base/case,pin)
    raw=base/case/'observations'/f"step_{receipt['terminal_step']:04d}"/'raw/image_0004'
    directory=tmp_path/'evaluation'; target=directory/'scf_000000'; target.mkdir(parents=True)
    for name in ('INPUT','KPT','STRU','OUT.ABACUS/running_scf.log'):
        destination=target/name; destination.parent.mkdir(exist_ok=True)
        shutil.copyfile(raw/name,destination)
    row=obs['raw_image_evaluations'][4]
    call=dict(input_sha256=copy.deepcopy(row['input_sha256']),raw_log_sha256=row['raw_log_sha256'],elapsed_seconds=1.)
    if fault=='input': call['input_sha256']['INPUT']='0'*64
    if fault=='zero_time': call['elapsed_seconds']=0.
    if fault=='nan_time': call['elapsed_seconds']=float('nan')
    if fault=='bool_time': call['elapsed_seconds']=True
    (target/'call_audit.json').write_text(json.dumps(call))
    real_sha=study.sha256
    monkeypatch.setattr(study,'sha256',lambda p: study.CONTRACT[p.name] if p.name in study.CONTRACT else real_sha(p))
    atoms=images[4]
    if fault in ('energy','force','stress'):
        raw_results=study.audited_results(target)
        if fault=='energy': raw_results['energy']+=1.
        if fault=='force': raw_results['forces'][0,0]+=.1
        if fault=='stress': raw_results['stress'][0]+=.1
        atoms.calc=SinglePointCalculator(atoms,**raw_results)
    if fault=='extra_SCF': (directory/'scf_000001').mkdir()
    if fault=='none':
        result=study.result_validator(0,atoms,directory)
        assert result['genuine_MPI_ranks']==32 and result['raw_SCFs']==1
        assert result['raw_log_sha256']==row['raw_log_sha256']
    else:
        with pytest.raises(ValueError): study.result_validator(0,atoms,directory)
