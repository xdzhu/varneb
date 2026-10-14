"""Finite canary preparation/orchestration; synthetic records are not DFT."""
from copy import deepcopy
import json
from pathlib import Path
import shutil

import numpy as np
import pytest

from scripts import hfo2_G3_conditional_canary as canary

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'benchmarks/hfo2_channels/20261008'


@pytest.fixture
def data():
    manifest=json.loads((BASE/'G3_curvature_E150_20261014/manifest.json').read_bytes())
    # Borrow the actual chart/shape only. This constructed report is explicitly
    # synthetic, not a native-completion artifact or a selected material plane.
    report=deepcopy(json.loads((BASE/'G3_first_amplitude_E155_20261014/T_native.json').read_bytes()))
    basis=np.asarray(manifest['centers']['T']['directions'])
    report.update(status='complete_local_curvature_diagnostic',completed_records=148,
        raw_SCFs_reparsed=148,native_record_sources=[{'explicit_mock_not_native':True} for _ in range(148)])
    internal=np.diag(np.r_[-2.,-.4,np.linspace(1.,3.,34)])
    reduced=np.zeros((37,37));reduced[:36,:36]=internal;reduced[-1,-1]=-7.
    report['blocks']=[dict(step_A=h,complete_pairs=37,full_stencil_at_this_step=True,
        raw_projected_eV_A2=reduced.tolist(),full_hessian_action_eV_A2=(basis@reduced).tolist(),
        fixed_strain_internal_eigenvalues_eV_A2=np.linalg.eigvalsh(internal).tolist()) for h in (.01,.02)]
    g=basis[:,3]*.1
    report['raw_center_gradient_eV_A']=g.tolist()
    return manifest,report


def test_actual_partial_prefix_refuses_before_selection_or_output():
    m=json.loads((BASE/'G3_curvature_E150_20261014/manifest.json').read_bytes())
    r=json.loads((BASE/'G3_first_amplitude_E155_20261014/T_native.json').read_bytes())
    with pytest.raises(ValueError,match='complete original two-amplitude'):
        canary.prepare(m,'T',r,'actual-partial-report-not-complete')


def test_mock_preparation_separates_external_control_and_does_not_alias_settings(data):
    m,r=data; c=canary.prepare(m,'T',r,'explicit_mock_not_native')
    assert c['raw_complement_minima_eV_A2']==pytest.approx([1.,1.])
    assert c['settings']['maximum_attempts_per_center']==24
    assert c['settings']['modeling_nodes_used']==1 and c['settings']['independent_test_points_used']==0
    assert not c['full_internal_index_or_TS_certified']
    c['settings']['arms'][0]='altered'
    assert canary.SETTINGS['arms']==['identity','curvature_seed']


def test_preserving_channel_uses_its_own_source_and_chart_without_new_data():
    m=json.loads((BASE/'G3_curvature_E150_20261014/manifest.json').read_bytes())
    r=json.loads((BASE/'G3_first_amplitude_E155_20261014/preserving_native.json').read_bytes())
    basis=np.asarray(m['centers']['preserving']['directions'])
    matrix=np.diag(np.r_[-2.,-.4,np.linspace(1.,3.,34),-7.])
    r.update(status='complete_local_curvature_diagnostic',completed_records=148,raw_SCFs_reparsed=148,
             native_record_sources=[{'explicit_mock_not_native':True} for _ in range(148)])
    r['blocks']=[dict(step_A=h,complete_pairs=37,full_stencil_at_this_step=True,
        raw_projected_eV_A2=matrix.tolist(),full_hessian_action_eV_A2=(basis@matrix).tolist(),
        fixed_strain_internal_eigenvalues_eV_A2=np.linalg.eigvalsh(matrix[:36,:36]).tolist()) for h in (.01,.02)]
    c=canary.prepare(m,'preserving',r,'explicit_mock_not_native')
    assert c['channel']=='preserving'
    canary.validate_config(c,m,'preserving',r,'explicit_mock_not_native')
    with pytest.raises(ValueError):canary.prepare(m,'T',r,'explicit_mock_not_native')


@pytest.mark.parametrize('fault',['missing_pair','raw_action','wrong_center','failure','one_negative','nonfinite'])
def test_mock_report_faults_do_not_generate_a_canary(data,fault):
    m,r=data
    if fault=='missing_pair':r['blocks'][1]['complete_pairs']=36
    elif fault=='raw_action':r['blocks'][0]['full_hessian_action_eV_A2'][0][0]+=.1
    elif fault=='wrong_center':r['center_source']={}
    elif fault=='failure':r['failures_present']=['explicit injected failure']
    elif fault=='nonfinite':r['blocks'][0]['raw_projected_eV_A2'][0][0]=np.nan
    else:
        basis=np.asarray(m['centers']['T']['directions'])
        for b in r['blocks']:
            a=np.asarray(b['raw_projected_eV_A2']);a[1,1]=.4
            b['raw_projected_eV_A2']=a.tolist();b['full_hessian_action_eV_A2']=(basis@a).tolist()
            b['fixed_strain_internal_eigenvalues_eV_A2']=np.linalg.eigvalsh(a[:36,:36]).tolist()
    with pytest.raises(ValueError):canary.prepare(m,'T',r,'explicit_mock_not_native')


@pytest.mark.parametrize('fault',['setting','bool_schema','code','axis_span','raw_minimum'])
def test_configuration_drift_is_refused(data,fault):
    m,r=data;c=canary.prepare(m,'T',r,'explicit_mock_not_native')
    if fault=='setting':c['settings']['maximum_attempts_per_arm']=13
    elif fault=='bool_schema':c['schema_version']=True
    elif fault=='code':c['normalized_runtime_source_sha256'][canary.CODE[0]]='0'*64
    elif fault=='raw_minimum':c['raw_complement_minima_eV_A2'][0]+=.1
    else:c['internal_axes_36_by_2']=np.eye(36)[:,2:4].tolist()
    with pytest.raises(ValueError):canary.validate_config(c,m,'T',r,'explicit_mock_not_native')


def test_prepared_axis_orientation_is_not_regenerated_by_sign_convention(data):
    m,r=data;c=canary.prepare(m,'T',r,'explicit_mock_not_native')
    axes=np.asarray(c['internal_axes_36_by_2']);c['internal_axes_36_by_2']=(-axes).tolist()
    canary.validate_config(c,m,'T',r,'explicit_mock_not_native')
    np.testing.assert_array_equal(c['internal_axes_36_by_2'],-axes)


def test_internal_pullback_reuses_only_exact_origin_and_preserves_raw_force_gradient(data):
    m,r=data;b=np.asarray(m['centers']['T']['directions'])[:,:36];g=np.asarray(r['raw_center_gradient_eV_A'])
    class Delegate:
        fixed_coordinate_covectors=np.eye(40)[:,-1:]
        def __init__(self):self.calls=[]
        def __call__(self,x):self.calls.append(x.copy());return 2.+g@x,g.copy()
    delegate=Delegate();wrapped=canary.InternalEvaluator(delegate,b,2.,g);z=np.zeros(36)
    e,gi=wrapped(z);gi[0]=999.
    assert not delegate.calls and wrapped.center_cache_hits==1
    z[2]=1e-14; e,gi=wrapped(z)
    assert len(delegate.calls)==1 and wrapped.center_cache_hits==1
    np.testing.assert_array_equal(delegate.calls[0],b@z)
    np.testing.assert_allclose(gi,b.T@g)
    assert wrapped.fixed_coordinate_covectors.shape==(36,0)
    np.testing.assert_array_equal(wrapped.source_fixed_coordinate_covectors,delegate.fixed_coordinate_covectors)
    assert wrapped.callback_requests==2 and not wrapped.gradient.flags.writeable
    with pytest.raises(ValueError):wrapped(np.zeros(36,dtype=complex))


@pytest.mark.parametrize('fault',['lost_control','mixed_basis','complex_gradient'])
def test_embedding_never_silently_drops_a_real_control(data,fault):
    m,r=data;b=np.asarray(m['centers']['T']['directions'])[:,:36].copy();g=np.asarray(r['raw_center_gradient_eV_A'])
    class Delegate:
        fixed_coordinate_covectors=np.eye(40)[:,-1:]
    delegate=Delegate()
    if fault=='lost_control':delegate.fixed_coordinate_covectors=np.zeros((40,0))
    elif fault=='mixed_basis':b[-1,0]=1e-13
    else:g=g.astype(complex)
    with pytest.raises(ValueError):canary.InternalEvaluator(delegate,b,0.,g)


def runtime(tmp_path,monkeypatch,data,mode='success'):
    m,r=data
    source=tmp_path/'source';source.mkdir()
    for name in canary.CODE:
        target=source/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/name,target)
    monkeypatch.setattr(canary.study,'ROOT',source);monkeypatch.setattr(canary.study,'RUNS',tmp_path)
    chart,_,g,matrices=canary.measured_pair(m,'T',r)
    monkeypatch.setattr(canary.study,'native_preflight',lambda *a:(chart,(),0.,g))
    monkeypatch.setattr(canary.study,'make_factory',lambda **a:object())
    for k,v in dict(RUN_DFT='1',SLURM_JOB_PARTITION='hfacnormal01',SLURM_JOB_NUM_NODES='1',
        SLURM_NTASKS='1',SLURM_CPUS_PER_TASK='32',SLURM_JOB_ID='TEST-NOT-SLURM',
        OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1').items():monkeypatch.setenv(k,v)
    created=[]
    class AnalyticDelegate:
        """File-shape mock only: no ASE/backend/native validator is called."""
        def __init__(self,chart,factory,**kwargs):
            self.fixed_coordinate_covectors=chart.fixed_coordinate_covectors
            self.workdir=kwargs['workdir'];self.workdir.mkdir()
            self.cap=kwargs['max_evaluation_attempts'];self.n_evaluation_attempts=self.n_new_evaluations=0
            self.cache={};created.append(self)
        def __call__(self,x):
            if x.tobytes() in self.cache:return self.cache[x.tobytes()]
            if self.n_evaluation_attempts>=self.cap:raise RuntimeError('mode evaluation attempt budget exhausted; calculator was not called')
            i=self.n_evaluation_attempts;self.n_evaluation_attempts+=1
            directory=self.workdir/f'eval-{i:06d}-explicit-mock';directory.mkdir()
            if mode=='native_failure':
                (directory/'failure.json').write_text('{"explicit_mock_not_native":true}')
                raise RuntimeError('explicit injected backend failure, not a budget condition')
            value=float(.5*x@x+g@x);gradient=x+g
            (directory/'result.json').write_text(json.dumps({'validation':{'elapsed_seconds':1.,'explicit_mock_not_native':True}}))
            self.n_new_evaluations+=1;self.cache[x.tobytes()]=(value,gradient)
            return value,gradient
    monkeypatch.setattr(canary,'CalculatorModeEvaluator',AnalyticDelegate)
    config=canary.prepare(m,'T',r,'explicit_mock_not_native')
    return m,r,config,created


def test_analytic_runtime_integration_preserves_two_fresh_arms_and_center_cost(tmp_path,monkeypatch,data):
    m,r,c,created=runtime(tmp_path,monkeypatch,data)
    out=canary.run(c,m,'T',r,'explicit_mock_not_native',tmp_path/'T')
    assert len(created)==2 and out['all_arms_gradient_converged']
    assert out['new_center_SCFs']==0 and not out['same_branch_or_material_speedup_certified']
    assert out['new_probe_attempts']<=24
    for row in out['arms']:
        assert row['center_cache_hits']>=1 and row['new_completed_SCFs']<=12
        assert row['callback_requests']==row['point']['callback_evaluations']
        assert row['point']['orthogonal_gradient_norm_eV_A']<=.03


def test_both_budget_outcomes_are_retained_without_marking_a_surface_complete(tmp_path,monkeypatch,data):
    m,r,c,created=runtime(tmp_path,monkeypatch,data)
    def exhaust(plane,q,evaluate,**kwargs):
        evaluate(np.zeros(36))
        for i in range(1,14):
            z=np.zeros(36);z[2]=i*1e-4
            assert kwargs['trial_validator'](z);evaluate(z)
    monkeypatch.setattr(canary,'relax_orthogonal_at_q',exhaust)
    out=canary.run(c,m,'T',r,'explicit_mock_not_native',tmp_path/'T')
    assert len(created)==2 and out['new_probe_attempts']==24
    assert not out['all_arms_gradient_converged']
    assert all(row['point'] is None and row['new_completed_SCFs']==12 for row in out['arms'])


def test_real_backend_error_category_is_not_treated_as_budget_or_retried(tmp_path,monkeypatch,data):
    m,r,c,created=runtime(tmp_path,monkeypatch,data,'native_failure')
    with pytest.raises(RuntimeError,match='explicit injected backend failure'):
        canary.run(c,m,'T',r,'explicit_mock_not_native',tmp_path/'T')
    assert len(created)==1 and len(list((tmp_path/'T').glob('identity/eval-*/failure.json')))==1
    assert not (tmp_path/'T/curvature_seed').exists()


def test_login_never_invokes_native_preflight(tmp_path,monkeypatch,data):
    m,r=data;c=canary.prepare(m,'T',r,'explicit_mock_not_native');calls=[]
    monkeypatch.delenv('RUN_DFT',raising=False)
    monkeypatch.setattr(canary.study,'native_preflight',lambda *a:calls.append(True))
    with pytest.raises(ValueError,match='allocation'):canary.run(c,m,'T',r,'explicit_mock_not_native',tmp_path/'T')
    assert not calls and not (tmp_path/'T').exists()


def test_native_preflight_route_is_read_only_and_checks_the_origin(tmp_path,monkeypatch,data):
    m,r,c,created=runtime(tmp_path,monkeypatch,data)
    monkeypatch.delenv('RUN_DFT',raising=False)
    chart,basis,energy,gradient,matrices=canary.native_canary_preflight(c,m,'T',r,'explicit_mock_not_native')
    assert not created and not (tmp_path/'T').exists()
    assert basis.shape==(40,36) and len(matrices)==2
    monkeypatch.setattr(canary,'geometry_guard',lambda *a:False)
    with pytest.raises(ValueError,match='center violates'):
        canary.native_canary_preflight(c,m,'T',r,'explicit_mock_not_native')
    assert not created


def test_existing_configuration_refuses_before_reading_material_inputs(tmp_path,monkeypatch):
    output=tmp_path/'configuration.json';output.write_text('original')
    monkeypatch.setattr('sys.argv',['canary','prepare','--parent','unused','--native-audit','unused',
        '--native-audit-sha256','unused','--channel','T','--config',str(output)])
    calls=[];monkeypatch.setattr(canary,'pinned',lambda *a:calls.append(True))
    with pytest.raises(FileExistsError):canary.main()
    assert not calls and output.read_text()=='original'


def test_scheduler_wrapper_keeps_native_resources_without_dispatch_or_alternate_launcher():
    text=(ROOT/'cluster/hf_hfo2_G3_conditional_canary.slurm').read_text()
    for line in ('#SBATCH -p hfacnormal01','#SBATCH -N 1','#SBATCH --ntasks=1',
                 '#SBATCH --cpus-per-task=32','#SBATCH --time=02:00:00'):
        assert line in text
    assert 'sbatch ' not in text and 'srun ' not in text
    assert 'PYTHONDONTWRITEBYTECODE=1' in text
    assert canary.study.COMMAND.startswith('mpirun -np 32 ')
