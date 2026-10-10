"""Single-use, case-specific G2 wave; no monitor, retry or parameter changes."""
import argparse
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile


RUNS = Path('/public/home/iai806/abacus/agent-runs')
ROOT = RUNS/'20261010-varneb-clamped-G2-E058-r1'
SOURCE = RUNS/'20261010-varneb-nested-controls-E055-r1/source'
ARCHIVE = RUNS/'varneb-E055-clean-86d2637.tar'
PREPARED = RUNS/'20261010-varneb-clamped-G2-remaining-E056-r1'
PARENTS = ('28661019', '28661020')
CHANNELS = ('PO_to_T', 'PO_flip_T_pattern_reversing')
DEPENDENCY = 'afterok:' + ':'.join(PARENTS)
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
    """Only the registered parents may be live in this study, not other users' jobs."""
    rows = [line.strip().split('|') for line in acct_text.splitlines() if line.strip()]
    if len(rows) != 2 or any(len(r) != 6 for r in rows):
        raise ValueError('missing or malformed parent accounting; reconcile, do not submit')
    acct = {r[0]: r for r in rows}
    if set(acct) != set(PARENTS):
        raise ValueError('unexpected parent handles')
    for job, state, exit_code, cpus, partition, user in rows:
        if (state not in ('RUNNING', 'COMPLETING', 'COMPLETED') or exit_code != '0:0'
                or cpus != '32' or partition != 'hfacnormal01' or user != 'iai806'):
            raise ValueError('parent is not a healthy registered 32-CPU allocation')
    study = []
    for line in queue_text.splitlines():
        fields = [f.strip() for f in line.split('|')]
        if len(fields) != 5:
            raise ValueError('malformed active queue row')
        if fields[1].startswith('hfo2'):
            study.append(fields)
    if len(study) > 2 or any(r[0] not in PARENTS or r[3:] != ['hfacnormal01', '32'] for r in study):
        raise ValueError('other study handles exist; do not duplicate or exceed concurrency')
    return dict(parent_rows=rows, active_study_rows=study, unrelated_jobs_untouched=True)


def scheduler_snapshot():
    acct = run(['sacct', '-X', '-n', '-j', ','.join(PARENTS), '-P',
                '--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,User'])
    queue = run(['squeue', '-h', '-u', 'iai806', '-o', '%i|%100j|%T|%P|%C'])
    result = reconcile(acct, queue)
    result.update(check_CST=datetime.now().isoformat(), accounting_raw=acct, queue_raw=queue)
    return result


def preflight():
    if ROOT.exists():
        raise FileExistsError('single-use namespace exists; read its receipt, never rerun blindly')
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
    shutil.copyfile(__file__, ROOT/'executed_queue_next_pair.py')
    import numpy as np
    from ase.io import read
    import examples.hfo2_fixed_input_factory as transport
    from vcneb import validate_path_geometry, validate_periodic_path_lift
    receipt = load(PREPARED/'preparation_receipt.json')
    expected = {r['channel']: r for r in receipt['prepared_channels'] if r['condition'] == 'strain_0000'}
    reports = []
    # The calculator may read actual raw outputs but cannot launch any executable here.
    real_run = subprocess.run
    def forbidden(*args, **kwargs):
        raise RuntimeError('preflight prohibits DFT/external executable launches')
    subprocess.run = forbidden
    try:
        for channel in CHANNELS:
            seed = PREPARED/'seeds/strain_0000'/channel
            m = load(seed/'manifest.json')
            if (sha(seed/'manifest.json') != expected[channel]['manifest_sha256']
                    or sha(seed/'vcneb_preflight_HF.json') != expected[channel]['geometry_preflight_sha256']
                    or any(sha(seed/n) != h for n,h in m['files_sha256'].items())
                    or m['strain'] != 0 or m['physical_contract_sha256'] != transport.CONTRACT
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
            parameters = load(seed/'factory_parameters.json')
            factory = transport.make_clamped_seed_cached_factory(parameters=parameters, command=COMMAND)
            caches = []
            for i, a in enumerate(images):
                directory = ROOT/'cache_preflight'/channel/f'image_{i:04d}'
                a.calc = factory(i, a, directory)
                if i in (0,8):
                    e, f, s = a.get_potential_energy(), a.get_forces(), a.get_stress()
                    if not np.isfinite(e) or f.shape != (12,3) or s.shape != (6,) or a.calc.next_call != 0:
                        raise ValueError('cached endpoint full EFS invalid')
                    caches.append(dict(image_index=i, energy_eV_cell=e, new_DFT_calls=0,
                        audit_sha256=sha(directory/'seed_cache_audit.json')))
                elif a.calc.results or a.calc.next_call != 0:
                    raise ValueError('interior must be fresh')
            reports.append(dict(channel=channel, seed_root=str(seed), manifest_sha256=sha(seed/'manifest.json'),
                seed_files_sha256=m['files_sha256'], endpoint_caches=caches, internal_images_fresh=7,
                physical_bytes_native_32MPI_raw_EFS_rechecked=True, geometry_checked_images=9))
    finally:
        subprocess.run = real_run
    result = dict(status='preflight_passed_NOT_submitted', check_CST=datetime.now().isoformat(),
        dependency=DEPENDENCY, scheduler=scheduler, registered_channels=list(CHANNELS), reports=reports,
        source_root=str(SOURCE), tested_source_commit='86d2637ebcd6c502df43b91ca03977fd50120716',
        archive_sha256=ARCHIVE_SHA, pilot_sha256=PILOT_SHA, code_files_byte_checked=len(code_hashes),
        source_code_sha256=code_hashes, queue_helper_sha256=sha(__file__), ABACUS_binary_sha256=BINARY_SHA,
        new_DFT_calls=0, existing_jobs_mutated=False, physical_parameters_changed=False,
        held_out_condition_generated_or_read=False, recurrence_or_retry_created=False,
        job_completion_is_not_NEB_convergence=True)
    with (ROOT/'preflight.json').open('x') as f:
        json.dump(result, f, indent=2); f.write('\n')
    return result


def sbatch_argv(channel):
    if channel not in CHANNELS:
        raise ValueError('only the two remaining zero-strain channels are authorized in this wave')
    return ['sbatch', '--parsable', '--dependency='+DEPENDENCY, '--job-name=hfo2-G2-E058-'+channel,
        '--export=ALL,RUN_DFT=1,SOURCE_ROOT='+SOURCE.as_posix()+',SEED_ROOT='+ (PREPARED/'seeds/strain_0000'/channel).as_posix()
        +',WORKDIR='+ (ROOT/channel/'band').as_posix(), '--output='+ (ROOT/(channel+'.slurm.out')).as_posix(),
        '--error='+ (ROOT/(channel+'.slurm.err')).as_posix(), (SOURCE/PILOT).as_posix()]


def submit_once():
    receipt = load(ROOT/'preflight.json')
    if (receipt['status'] != 'preflight_passed_NOT_submitted' or sha(__file__) != receipt['queue_helper_sha256']
            or sha(SOURCE/PILOT) != PILOT_SHA or receipt['registered_channels'] != list(CHANNELS)):
        raise ValueError('reviewed submission source/evidence changed')
    for r in receipt['reports']:
        seed = Path(r['seed_root'])
        if sha(seed/'manifest.json') != r['manifest_sha256'] or any(sha(seed/n) != h for n,h in r['seed_files_sha256'].items()):
            raise ValueError('seed changed since preflight')
        if (ROOT/r['channel']).exists():
            raise FileExistsError('existing production workdir')
    scheduler = scheduler_snapshot()
    journal_path = ROOT/'submission_journal.jsonl'
    jobs = []
    # Exclusive journal creation also protects a partial/ambiguous first submission.
    with journal_path.open('x') as journal:
        def record(item):
            journal.write(json.dumps(item)+'\n'); journal.flush(); os.fsync(journal.fileno())
        record(dict(event='reviewed_intent', check_CST=datetime.now().isoformat(), dependency=DEPENDENCY,
                    channels=list(CHANNELS), preflight_sha256=sha(ROOT/'preflight.json'), scheduler=scheduler))
        for channel in CHANNELS:
            argv = sbatch_argv(channel)
            record(dict(event='submission_attempt', channel=channel, argv=argv))
            try:
                stdout = run(argv).strip()
                if not stdout.isdecimal():
                    raise ValueError('ambiguous sbatch reply: ' + stdout)
            except Exception as exc:
                record(dict(event='ambiguous_or_failed_reply_RECONCILE_DO_NOT_RETRY', channel=channel,
                            exception_type=type(exc).__name__, message=str(exc)))
                raise
            record(dict(event='accepted_handle', channel=channel, job_id=stdout, dependency=DEPENDENCY))
            jobs.append(dict(channel=channel, job_id=stdout))
    result = dict(status='two_registered_starts_submitted_with_dependency', check_CST=datetime.now().isoformat(),
        jobs=jobs, dependency=DEPENDENCY, preflight_sha256=sha(ROOT/'preflight.json'),
        submission_journal_sha256=sha(journal_path), parent_jobs_not_modified=True,
        max_study_simultaneous_allocations=2, total_registered_G2_matrix_chains=8,
        G2_unique_chains_with_handles=4, G2_training_starts_remaining_without_handles=4,
        physical_parameters_changed=False, held_out_condition_generated_or_read=False,
        segment_steps=10, segment_walltime_hours=4, job_completion_is_not_NEB_convergence=True)
    with (ROOT/'submission_receipt.json').open('x') as f:
        json.dump(result, f, indent=2); f.write('\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--preflight', action='store_true')
    group.add_argument('--submit', action='store_true')
    args = parser.parse_args()
    result = preflight() if args.preflight else submit_once()
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_code_sha256','reports')}))
