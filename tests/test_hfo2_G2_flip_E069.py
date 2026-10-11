import importlib.util
import json
from pathlib import Path
import pytest

from scripts.audit_hfo2_clamped_terminal import validate_terminal
from scripts.audit_hfo2_static_replica import sha256
from scripts.export_hfo2_clamped_observation import REGISTERED_FACTORIES,FLIP_E069_SCRIPT_SHA256

REPO = Path(__file__).resolve().parents[1]
CASE = REPO/'benchmarks/hfo2_channels/20261008/clamped_flip_continue_E069_20261011'


def load(path,name):
    spec = importlib.util.spec_from_file_location(name,path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def cap():
    return dict(status='max_steps_reached',converged=False,termination='step_limit',
        final_max_generalized_force_eV_per_A=.197683006020693,climbing_image_active_final=False)


def test_default_ordinary_gate_not_weakened_by_explicit_stepcap_mode():
    row = '28722810|COMPLETED|0:0|32|hfacnormal01|iai806'
    with pytest.raises(ValueError,match='ordinary'):
        validate_terminal('28722810',row,cap(),False)
    validate_terminal('28722810',row,cap(),False,outcome='step-cap')
    with pytest.raises(ValueError,match='outcome'):
        validate_terminal('28722810',row,cap(),False,outcome='guess')


@pytest.mark.parametrize('change',[dict(status='failed'),dict(converged=True),dict(termination='error'),
    dict(final_max_generalized_force_eV_per_A=.10),dict(final_max_generalized_force_eV_per_A=-1),
    dict(final_max_generalized_force_eV_per_A=True),dict(final_max_generalized_force_eV_per_A=float('nan')),
    dict(final_max_generalized_force_eV_per_A=float('inf')),dict(climbing_image_active_final=True)])
def test_no_failed_cap_nonfinite_or_converged_seed(change):
    with pytest.raises(ValueError,match='step cap'):
        validate_terminal('28722810','28722810|COMPLETED|0:0|32|hfacnormal01|iai806',
            {**cap(),**change},False,outcome='step-cap')


def test_successful_cap_still_refuses_failure_record_or_live_accounting():
    with pytest.raises(ValueError,match='step cap'):
        validate_terminal('28722810','28722810|COMPLETED|0:0|32|hfacnormal01|iai806',cap(),True,outcome='step-cap')
    with pytest.raises(ValueError,match='actual completed'):
        validate_terminal('28722810','28722810|RUNNING|0:0|32|hfacnormal01|iai806',cap(),False,outcome='step-cap')


def test_exact_E069_recipe_not_a_physical_input_change():
    path = REPO/'cluster/hf_hfo2_clamped_flip_continue_E069_20261011.slurm'
    assert sha256(path) == FLIP_E069_SCRIPT_SHA256
    assert REGISTERED_FACTORIES[FLIP_E069_SCRIPT_SHA256].endswith(':make_clamped_resume_cached_factory')
    text = path.read_text()
    for token in ('#SBATCH -p hfacnormal01','--cpus-per-task=32','--time=12:00:00',
                  'source_job_id":"28722810"','--steps 40','--fmax .10','--no-climb','--pressure-gpa 0',
                  'mpirun -np 32','--maxstep .02 --k .2','source-fixed','I_MPI_HYDRA_BOOTSTRAP=fork'):
        assert token in text
    assert not any(token in text for token in ('ecutwfc','ecutrho','srun ','SYMPREC','--optimizer LBFGS'))


@pytest.mark.parametrize('extra',[None,'unknown','second_active','wrong_cpu','wrong_partition','malformed'])
def test_capacity_ignores_other_project_but_never_another_study_slot(extra):
    module = load(CASE/'dispatch_once.py','E069_queue')
    lines = ['28692775|hfo2-G2-E058-PO_to_T|RUNNING|hfacnormal01|32',
             '999|pc-unrelated|RUNNING|another|64']
    if extra == 'unknown': lines += ['998|hfo2-other|PENDING|hfacnormal01|32']
    if extra == 'second_active': lines += ['28692776|hfo2-G2|RUNNING|hfacnormal01|32']
    if extra == 'wrong_cpu': lines[0] = lines[0][:-2]+'1'
    if extra == 'wrong_partition': lines[0] = lines[0].replace('hfacnormal01','other')
    if extra == 'malformed': lines += ['bad']
    if extra:
        with pytest.raises(ValueError): module.spare_slot('\n'.join(lines))
    else:
        assert module.spare_slot('\n'.join(lines)) == ['28692775']


@pytest.fixture
def reviewed(tmp_path,monkeypatch):
    module = load(CASE/'dispatch_once.py','E069_dispatch')
    parser = load(REPO/'benchmarks/hfo2_channels/20261008/clamped_flip_continue_E067_20261011/submit_held_once.py','E067_pure_parser')
    module.ROOT = tmp_path
    module.SCRIPT = tmp_path/'script.slurm'
    module.SCRIPT.write_text('# reviewed\n')
    monkeypatch.setattr(module,'verify_audit',lambda digest: dict(seed_manifest_sha256='a'*64))
    monkeypatch.setattr(module,'rules',lambda: parser)
    state = dict(writes=[],ambiguous=None,accepted_bad=False)
    def invoke(argv):
        if argv[0] == 'sbatch' or argv[:2] == ['scontrol','release']:
            state['writes'].append(argv)
            if state['ambiguous'] == argv[0]: raise TimeoutError('ambiguous external response')
            return '28723000' if argv[0] == 'sbatch' else ''
        if argv[0] == 'sacct': return '28722810|COMPLETED|0:0|32|hfacnormal01|iai806'
        if argv[0] == 'squeue':
            return '28692775|hfo2-G2-E058-PO_to_T|RUNNING|hfacnormal01|32'
        if argv[:3] == ['scontrol','show','job']:
            job = argv[-1]
            channel = 'PO_flip_T_pattern_reversing'
            if job == '28692776':
                row = dict(JobId=job,JobName='hfo2-G2-E058-'+channel,UserId='iai806(1)',JobState='PENDING',
                    Reason='JobHeldUser',Priority='0',Partition='hfacnormal01',NumCPUs='32',Dependency='(null)',
                    Command=parser.PILOT,StdOut=parser.WAIT_ROOT+'/'+channel+'.slurm.out',StdErr=parser.WAIT_ROOT+'/'+channel+'.slurm.err')
            else:
                row = dict(JobId=job,JobName='hfo2-G2-flip-E069',UserId='iai806(1)',JobState='PENDING',
                    Reason='JobHeldUser',Priority='0',Partition='hfacnormal01',NumCPUs='32',Dependency='(null)',
                    Command=str(module.SCRIPT),StdOut=str(tmp_path/'flip.slurm.out'),StdErr=str(tmp_path/'flip.slurm.err'))
                if state['accepted_bad']: row['NumCPUs']='1'
            return ' '.join(f'{k}={v}' for k,v in row.items())
        if argv[0] == 'bash': return ''
        raise AssertionError(argv)
    monkeypatch.setattr(module,'run',invoke)
    return module,state


def test_one_same_chain_dispatch_and_no_dependency_writes(reviewed):
    module,state = reviewed
    result = module.dispatch_once('b'*64)
    assert [v[0] for v in state['writes']] == ['sbatch','scontrol']
    assert state['writes'][1] == ['scontrol','release','28723000']
    assert result['total_submissions'] == result['release_writes'] == 1
    assert result['dependency_writes'] == result['new_independent_chains'] == 0
    assert result['source_step'] == 20 and result['max_new_updates'] == 40
    assert not result['execution_started_claimed']
    assert result['journal_sha256'] == module.sha(module.ROOT/'dispatch_journal.jsonl')
    with pytest.raises(FileExistsError): module.dispatch_once('b'*64)
    assert len(state['writes']) == 2


@pytest.mark.parametrize('failure',['sbatch','scontrol','bad_identity'])
def test_ambiguous_or_bad_accepted_state_never_repeats_submission_or_release(reviewed,failure):
    module,state = reviewed
    state['ambiguous'] = failure
    state['accepted_bad'] = failure == 'bad_identity'
    with pytest.raises((TimeoutError,ValueError)): module.dispatch_once('b'*64)
    writes = list(state['writes'])
    with pytest.raises(FileExistsError): module.dispatch_once('b'*64)
    assert writes == state['writes']
    assert len(writes) == (2 if failure == 'scontrol' else 1)
    assert not (module.ROOT/'dispatch_receipt.json').exists()
    assert 'reconcile_existing_intent_never_repeat_writes' in (module.ROOT/'dispatch_journal.jsonl').read_text()


def test_actual_native_terminal_and_exact_latest_seed_replay():
    import numpy as np
    from ase.io import read
    from ase.calculators.singlepoint import SinglePointCalculator
    from examples.hfo2_fixed_input_factory import CONTRACT,same_ordered_geometry,read_fixed_hfo2_stru
    from scripts.audit_hfo2_static_replica import audited_results
    from scripts.analyze_hfo2_clamped_residual import analyze
    from vcneb import VCNEB,clamped_plane_vcneb_boundary
    receipt = json.loads((CASE/'audit_receipt.json').read_text())
    assert sha256(CASE/'audit_receipt.json') == '5f078aca099f9004949eaee7dd28307adec9fdced6300b592243a7d1138802e7'
    validate_terminal('28722810',receipt['scheduler_row'],json.loads((CASE/'vcneb_summary.json').read_text()),False,outcome='step-cap')
    assert sha256(CASE/'vcneb_summary.json') == receipt['terminal_summary_sha256']
    assert sha256(CASE/'executed_auditor.py') == receipt['executed_auditor_sha256']
    assert sha256(CASE/'executed_exporter.py') == receipt['executed_exporter_sha256']
    assert receipt['ordinary_converged'] is False and receipt['terminal_step'] == 20
    assert receipt['fresh_interior_SCFs'] == len(receipt['calls']) == 140
    assert len(receipt['runtime_code_sha256']) == 334
    assert receipt['new_DFT_calls'] == receipt['scheduler_mutations'] == 0
    assert receipt['current_frame_exact_caches_checked'] == 9
    assert receipt['fresh_SCF_wall_seconds_sum'] == pytest.approx(12585.1218851171,abs=1e-9,rel=0)
    r = receipt['terminal_observation']
    directory = CASE/'observations/step_0020'
    assert sha256(directory/'evaluated_chain.traj') == r['evaluated_chain_sha256']
    images = read(directory/'evaluated_chain.traj',index=':')
    for i,image in enumerate(images):
        raw = directory/'raw'/f'image_{i:04d}'
        p = r['raw_image_evaluations'][i]
        assert all(p['input_sha256'][n] == h for n,h in CONTRACT.items())
        for n in ('INPUT','KPT','STRU'): assert sha256(raw/n) == p['input_sha256'][n]
        assert sha256(raw/'OUT.ABACUS/running_scf.log') == p['raw_log_sha256']
        assert same_ordered_geometry(image,read_fixed_hfo2_stru(raw/'STRU'))
        parsed = audited_results(raw)
        for key,expected in (('energy',p['energy_eV_cell']),('forces',p['forces_eV_A']),('stress',p['stress_ASE_voigt_eV_A3'])):
            np.testing.assert_allclose(parsed[key],expected,atol=1e-12,rtol=0)
        cache = json.loads((CASE/'cache_preflight'/f'image_{i:04d}'/'seed_cache_audit.json').read_text())
        assert cache['policy'] == 'identical_ordered_clamped_resume_hash_pinned'
        assert cache['new_DFT_calls'] == 0 and cache['input_sha256'] == p['input_sha256']
        assert cache['raw_log_sha256'] == p['raw_log_sha256']
        image.calc = SinglePointCalculator(image,**parsed)
    boundary = clamped_plane_vcneb_boundary(12,np.asarray(r['mechanical_boundary']['reference_cell_A']),allow_tilt=True)
    band = VCNEB(images,cell_scale=r['cell_scale_A'],k=.2,climb=False,pressure=0,**boundary.vcneb_kwargs(images))
    assert np.linalg.norm(band.get_forces(),axis=1).max() == pytest.approx(.197683006020693,abs=1e-10,rel=0)
    profile = (band.enthalpies-band.enthalpies[0])*1000/4
    np.testing.assert_allclose(profile,r['relative_enthalpy_meV_fu'],atol=1e-9,rtol=0)
    assert profile.max() == pytest.approx(180.008886896303,abs=1e-9,rel=0)
    assert analyze(directory)['ordinary_residual_pass'] is False
    seed = CASE/'seed'
    manifest = json.loads((seed/'manifest.json').read_text())
    assert manifest['source_job_id'] == '28722810' and manifest['source_step'] == 20
    assert sha256(seed/'manifest.json') == receipt['seed_manifest_sha256']
    assert all(sha256(seed/n) == h for n,h in manifest['files_sha256'].items())
    assert sha256(seed/'seed.traj') == r['evaluated_chain_sha256']


def test_actual_single_dispatch_and_independent_native_start():
    receipt = json.loads((CASE/'dispatch_receipt.json').read_text())
    journal = [json.loads(line) for line in (CASE/'dispatch_journal.jsonl').read_text().splitlines()]
    assert receipt['audit_receipt_sha256'] == sha256(CASE/'audit_receipt.json')
    assert receipt['executed_helper_sha256'] == sha256(CASE/'dispatch_once.py')
    assert receipt['journal_sha256'] == sha256(CASE/'dispatch_journal.jsonl')
    assert receipt['total_submissions'] == receipt['release_writes'] == 1
    assert receipt['dependency_writes'] == receipt['new_independent_chains'] == 0
    assert receipt['job_id'] == '28723655' and receipt['parent_job'] == '28722810'
    assert [r['event'] for r in journal] == ['one_held_submission_intent','one_accepted_handle','one_release_intent','one_existing_handle_released']
    assert all(r['job_id'] == '28723655' for r in journal[1:])
    assert receipt['max_new_updates'] == 40 and not receipt['execution_started_claimed']
    startup = json.loads((CASE/'startup_observation.json').read_text())
    assert startup['raw_accounting'] == '28723655|RUNNING|0:0|32|2026-10-11T10:41:39|hfacnormal01'
    assert 'DSIZE = 32' in startup['native_first_fresh_SCF_DSIZE_line']
    assert startup['study_running_allocations'] == 2 and startup['stderr_bytes'] == 0
    assert not startup['fresh_SCF_completion_claimed']
