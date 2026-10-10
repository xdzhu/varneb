"""Single-use four-start training queue; no monitor, retry or DFT preflight."""
import argparse
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tarfile


RUNS = Path('/public/home/iai806/abacus/agent-runs')
ROOT = RUNS/'20261010-varneb-clamped-G2-E060-r1'
SOURCE = RUNS/'20261010-varneb-nested-controls-E055-r1/source'
ARCHIVE = RUNS/'varneb-E055-clean-86d2637.tar'
PREPARED = RUNS/'20261010-varneb-clamped-G2-remaining-E056-r1'
KNOWN = ('28661019', '28661020', '28692775', '28692776')
PREDECESSORS = KNOWN[2:]
WAVES = (('PO_to_T', 'PO_to_M'),
         ('PO_flip_T_pattern_preserving', 'PO_flip_T_pattern_reversing'))
CHANNELS = tuple(c for pair in WAVES for c in pair)
CONDITION = 'strain_p0100'
PILOT = 'cluster/hf_hfo2_clamped_chain_pilot_20261010.slurm'
ARCHIVE_SHA = 'cbe6552a3fa4237fec2098c08c8056735a9a95df4e18691d09fed6b4f24a613c'
PILOT_SHA = '7bf764238c4cc857ca2b7e3e94642f0517e1f3b84f4e4c5b17833379ac1101ed'
PREPARATION_SHA = '2cf1b83e32d5c24dfe99a981d5821f608533b51f4a43933e46e840fcc8c9c826'
BINARY = '/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus'
BINARY_SHA = '88ea7f91c9f3091ef5fb4f410ab88267b14ecfe7f35c86dfaee9bf9319cace7a'
COMMAND = 'mpirun -np 32 ' + BINARY


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def run(argv):
    return subprocess.run(argv, check=True, capture_output=True, text=True, timeout=45).stdout


def reconcile(acct_text, queue_text):
    rows = [line.strip().split('|') for line in acct_text.splitlines() if line.strip()]
    if len(rows) != 4 or any(len(r) != 6 for r in rows) or {r[0] for r in rows} != set(KNOWN):
        raise ValueError('all four existing study handles require authoritative reconciliation')
    for job, state, exit_code, cpus, partition, user in rows:
        expected_cpus = '0' if state == 'PENDING' else '32'
        if (state not in ('PENDING', 'RUNNING', 'COMPLETING', 'COMPLETED')
                or exit_code != '0:0' or cpus != expected_cpus
                or partition != 'hfacnormal01' or user != 'iai806'):
            raise ValueError('existing study handle failed or its ownership/resources changed')
    study = []
    for line in queue_text.splitlines():
        if not line.strip():
            continue
        fields = [f.strip() for f in line.split('|')]
        if len(fields) != 5:
            raise ValueError('malformed active queue row')
        if fields[1].startswith('hfo2'):
            study.append(fields)
    if (len({r[0] for r in study}) != len(study)
            or any(r[0] not in KNOWN or r[2] not in ('PENDING', 'RUNNING', 'COMPLETING')
                   or r[3:] != ['hfacnormal01', '32'] for r in study)
            or sum(r[2] != 'PENDING' for r in study) > 2):
        raise ValueError('unknown/duplicate study handle or concurrency violation')
    return dict(existing_accounting_rows=rows, study_queue_rows=study,
                active_study_allocations=sum(r[2] != 'PENDING' for r in study),
                unrelated_jobs_untouched=True)


def validate_predecessor_controls(rows, controls):
    acct = {r[0]:r for r in rows}
    remaining = {j for j in KNOWN[:2] if acct[j][1] != 'COMPLETED'}
    if set(controls) != set(PREDECESSORS):
        raise ValueError('both zero-strain predecessor controls required')
    evidence = []
    for i,job in enumerate(PREDECESSORS):
        fields = dict(t.split('=',1) for t in controls[job].split() if '=' in t)
        if (fields.get('JobId') != job or fields.get('JobState') != acct[job][1]
                or fields.get('JobName') != 'hfo2-G2-E058-'+('PO_to_T','PO_flip_T_pattern_reversing')[i]
                or fields.get('Partition') != 'hfacnormal01' or fields.get('NumCPUs') != '32'
                or fields.get('NumTasks') != '1' or fields.get('CPUs/Task') != '32'
                or fields.get('NumNodes') not in ('1','1-1')
                or fields.get('Command') != (SOURCE/PILOT).as_posix()):
            raise ValueError('existing zero-strain source, resources or state changed')
        dep = fields.get('Dependency')
        if dep is None:
            raise ValueError('dependency must be observed, not assumed')
        tokens = [] if dep == '(null)' else dep.split(',')
        parsed = [re.fullmatch(r'afterok:(\d+)(?:\((?:unfulfilled|fulfilled)\))?',t) for t in tokens]
        if any(m is None for m in parsed):
            raise ValueError('unexpected dependency type or expression')
        handles = [m.group(1) for m in parsed]
        if len(set(handles)) != len(handles) or not set(handles).issubset(KNOWN[:2]):
            raise ValueError('unexpected dependency handle')
        if fields['JobState'] == 'PENDING':
            if not remaining.issubset(handles):
                raise ValueError('zero-strain predecessors can escape active parent dependency')
        elif remaining:
            raise ValueError('zero-strain successor active before registered parents completed')
        evidence.append(dict(job_id=job,actual_fields=fields,dependency_graph_verified=True))
    return evidence


def scheduler_snapshot():
    acct = run(['sacct', '-X', '-n', '-j', ','.join(KNOWN), '-P',
                '--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,User'])
    queue = run(['squeue', '-h', '-u', 'iai806', '-o', '%i|%100j|%T|%P|%C'])
    result = reconcile(acct, queue)
    controls = {j:run(['scontrol','show','job','-o',j]) for j in PREDECESSORS}
    result['zero_strain_predecessor_controls'] = validate_predecessor_controls(result['existing_accounting_rows'],controls)
    result['zero_strain_control_raw'] = controls
    result.update(check_CST=datetime.now().isoformat(), accounting_raw=acct, queue_raw=queue)
    return result


def preflight():
    if ROOT.exists():
        raise FileExistsError('single-use namespace exists; reconcile, do not rerun blindly')
    if sha(ARCHIVE) != ARCHIVE_SHA or sha(SOURCE/PILOT) != PILOT_SHA:
        raise ValueError('tested immutable archive/pilot changed')
    if sha(PREPARED/'preparation_receipt.json') != PREPARATION_SHA or sha(BINARY) != BINARY_SHA:
        raise ValueError('prepared evidence or ABACUS binary changed')
    code_hashes = {}
    with tarfile.open(ARCHIVE) as bundle:
        for entry in bundle.getmembers():
            if entry.isfile() and entry.name.endswith('.py') and entry.name.startswith(('vcneb/', 'scripts/', 'examples/')):
                raw = bundle.extractfile(entry).read()
                if raw != (SOURCE/entry.name).read_bytes():
                    raise ValueError('tested source changed: ' + entry.name)
                code_hashes[entry.name] = hashlib.sha256(raw).hexdigest()
    if len(code_hashes) != 338:
        raise ValueError('unexpected archive coverage')
    run(['bash', '-n', str(SOURCE/PILOT)])
    scheduler = scheduler_snapshot()
    ROOT.mkdir(exist_ok=False)
    shutil.copyfile(__file__, ROOT/'executed_queue_training_waves.py')
    import numpy as np
    from ase.io import read
    import examples.hfo2_fixed_input_factory as transport
    from vcneb import validate_path_geometry, validate_periodic_path_lift
    preparation = load(PREPARED/'preparation_receipt.json')
    expected = {r['channel']:r for r in preparation['prepared_channels'] if r['condition'] == CONDITION}
    if set(expected) != set(CHANNELS):
        raise ValueError('only the four registered +1% starts may be read')
    reports = []
    real_run = subprocess.run
    def forbidden(*args, **kwargs):
        raise RuntimeError('preflight prohibits DFT/external executable launches')
    subprocess.run = forbidden
    try:
        for channel in CHANNELS:
            seed = PREPARED/'seeds'/CONDITION/channel
            m = load(seed/'manifest.json')
            if (sha(seed/'manifest.json') != expected[channel]['manifest_sha256']
                    or sha(seed/'vcneb_preflight_HF.json') != expected[channel]['geometry_preflight_sha256']
                    or any(sha(seed/n) != h for n,h in m['files_sha256'].items())
                    or m['strain'] != .01 or m['physical_contract_sha256'] != transport.CONTRACT
                    or m['physical_inputs_changed'] or m['holdout_generated'] or m['climb']
                    or m['fmax_eV_A'] != .10 or m['channel'] != channel):
                raise ValueError('registered training seed changed')
            images = read(seed/'seed.traj', index=':')
            if len(images) != 9 or any(a.calc is not None for a in images):
                raise ValueError('nine calculator-free seed geometries required')
            validate_periodic_path_lift(images)
            validate_path_geometry(images, cell_scale=m['cell_scale_A'], minimum_distance=1.6, maximum_deformation=.25)
            for a in images:
                if (a.get_chemical_symbols() != ['Hf']*4+['O']*8 or not np.allclose(
                        a.cell.array[:2], np.asarray(m['reference_cell_A'])[:2], atol=1e-10, rtol=0)):
                    raise ValueError('ordered atoms or common substrate plane changed')
            factory = transport.make_clamped_seed_cached_factory(
                parameters=load(seed/'factory_parameters.json'), command=COMMAND)
            caches = []
            for i,a in enumerate(images):
                directory = ROOT/'cache_preflight'/channel/f'image_{i:04d}'
                a.calc = factory(i,a,directory)
                if i in (0,8):
                    e,f,s = a.get_potential_energy(),a.get_forces(),a.get_stress()
                    if not np.isfinite(e) or f.shape != (12,3) or s.shape != (6,) or a.calc.next_call != 0:
                        raise ValueError('cached endpoint full EFS invalid')
                    caches.append(dict(image_index=i,energy_eV_cell=e,new_DFT_calls=0,
                                       audit_sha256=sha(directory/'seed_cache_audit.json')))
                elif a.calc.results or a.calc.next_call != 0:
                    raise ValueError('interior must be fresh')
            reports.append(dict(channel=channel,condition=CONDITION,seed_root=str(seed),
                manifest_sha256=sha(seed/'manifest.json'),seed_files_sha256=m['files_sha256'],
                endpoint_caches=caches,internal_images_fresh=7,geometry_checked_images=9,
                physical_bytes_native_32MPI_raw_EFS_rechecked=True))
    finally:
        subprocess.run = real_run
    result = dict(status='four_training_starts_preflight_passed_NOT_submitted',check_CST=datetime.now().isoformat(),
        predecessor_handles=list(PREDECESSORS),waves=WAVES,condition=CONDITION,scheduler=scheduler,reports=reports,
        source_root=str(SOURCE),tested_source_commit='86d2637ebcd6c502df43b91ca03977fd50120716',
        archive_sha256=ARCHIVE_SHA,pilot_sha256=PILOT_SHA,code_files_byte_checked=len(code_hashes),
        source_code_sha256=code_hashes,queue_helper_sha256=sha(__file__),ABACUS_binary_sha256=BINARY_SHA,
        new_DFT_calls=0,existing_jobs_mutated=False,physical_parameters_changed=False,
        held_out_condition_generated_or_read=False,recurrence_or_retry_created=False)
    with (ROOT/'preflight.json').open('x') as f:
        json.dump(result,f,indent=2);f.write('\n')
    return result


def sbatch_argv(channel, predecessors):
    if channel not in CHANNELS:
        raise ValueError('only the four registered +1% training channels are authorized')
    if (not isinstance(predecessors,tuple) or len(predecessors) != 2
            or any(not isinstance(j,str) or not j.isdecimal() or int(j) <= 0 for j in predecessors)
            or len(set(predecessors)) != 2):
        raise ValueError('two distinct explicit predecessor handles required')
    dependency = 'afterok:' + ':'.join(predecessors)
    return ['sbatch','--parsable','--dependency='+dependency,'--job-name=hfo2-G2-E060-'+channel,
        '--export=ALL,RUN_DFT=1,SOURCE_ROOT='+SOURCE.as_posix()+',SEED_ROOT='+(PREPARED/'seeds'/CONDITION/channel).as_posix()
        +',WORKDIR='+(ROOT/channel/'band').as_posix(),'--output='+(ROOT/(channel+'.slurm.out')).as_posix(),
        '--error='+(ROOT/(channel+'.slurm.err')).as_posix(),(SOURCE/PILOT).as_posix()]


def submit_once():
    receipt = load(ROOT/'preflight.json')
    if (receipt['status'] != 'four_training_starts_preflight_passed_NOT_submitted'
            or sha(__file__) != receipt['queue_helper_sha256'] or sha(SOURCE/PILOT) != PILOT_SHA
            or receipt['waves'] != [list(p) for p in WAVES] or receipt['condition'] != CONDITION):
        raise ValueError('reviewed submission source/evidence changed')
    if {r['channel'] for r in receipt['reports']} != set(CHANNELS):
        raise ValueError('incomplete training preflight')
    for r in receipt['reports']:
        seed = Path(r['seed_root'])
        if sha(seed/'manifest.json') != r['manifest_sha256'] or any(sha(seed/n) != h for n,h in r['seed_files_sha256'].items()):
            raise ValueError('seed changed since preflight')
        if (ROOT/r['channel']).exists():
            raise FileExistsError('existing production workdir')
    scheduler = scheduler_snapshot()
    journal_path = ROOT/'submission_journal.jsonl'
    jobs = []
    with journal_path.open('x') as journal:
        def record(item):
            journal.write(json.dumps(item)+'\n');journal.flush();os.fsync(journal.fileno())
        record(dict(event='reviewed_intent',check_CST=datetime.now().isoformat(),waves=WAVES,
                    condition=CONDITION,preflight_sha256=sha(ROOT/'preflight.json'),scheduler=scheduler))
        predecessors = PREDECESSORS
        for pair in WAVES:
            next_handles = []
            dependency = 'afterok:' + ':'.join(predecessors)
            for channel in pair:
                argv = sbatch_argv(channel,predecessors)
                record(dict(event='submission_attempt',channel=channel,dependency=dependency,argv=argv))
                try:
                    reply = run(argv).strip()
                    if not reply.isdecimal() or int(reply) <= 0 or reply in KNOWN or any(j['job_id'] == reply for j in jobs):
                        raise ValueError('ambiguous or duplicate sbatch reply: '+reply)
                except Exception as exc:
                    record(dict(event='ambiguous_or_failed_reply_RECONCILE_DO_NOT_RETRY',channel=channel,
                                exception_type=type(exc).__name__,message=str(exc)))
                    raise
                record(dict(event='accepted_handle',channel=channel,job_id=reply,dependency=dependency))
                jobs.append(dict(channel=channel,condition=CONDITION,job_id=reply,dependency=dependency))
                next_handles.append(reply)
            predecessors = tuple(next_handles)
    result = dict(status='four_registered_training_starts_submitted_in_two_dependency_waves',
        check_CST=datetime.now().isoformat(),jobs=jobs,predecessor_handles=list(PREDECESSORS),
        preflight_sha256=sha(ROOT/'preflight.json'),submission_journal_sha256=sha(journal_path),
        max_study_simultaneous_allocations=2,total_registered_G2_matrix_chains=8,
        G2_unique_chains_with_handles=8,G2_training_starts_remaining_without_handles=0,
        segment_steps=10,segment_walltime_hours=4,parent_or_previous_jobs_modified=False,
        job_completion_is_not_NEB_convergence=True,physical_parameters_changed=False,
        held_out_condition_generated_or_read=False,recurrence_or_retry_created=False)
    with (ROOT/'submission_receipt.json').open('x') as f:
        json.dump(result,f,indent=2);f.write('\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--preflight',action='store_true')
    group.add_argument('--submit',action='store_true')
    args = parser.parse_args()
    result = preflight() if args.preflight else submit_once()
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_code_sha256','reports')}))
