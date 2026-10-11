"""Single reviewed rolling-slot handoff; never submit or alter physical inputs."""
from datetime import datetime
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time

RUNS = Path('/public/home/iai806/abacus/agent-runs')
ROOT = RUNS/'20261011-varneb-clamped-M-terminal-E068-r1'
DONE = '28722320'
LIVE = '28722810'
TARGET = '28692775'
OTHER_HELD = '28692776'
DOWNSTREAM = {'28709788','28709789','28709790','28709791'}
WAIT_ROOT = RUNS/'20261010-varneb-clamped-G2-E058-r1'
RULE_SOURCE = RUNS/'20261011-varneb-clamped-flip-E067-r1/submit_held_once.py'
RULE_SHA = 'e74ba9ba9090ec0a46cb56d046c4931f7889e4a9e4855a38fe59544b12708927'
PREFLIGHT_SHA = '069107a8cdaca498262396a3fe5be315a1c739ba5154a937fc742c099bb371fb'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_rules():
    if sha(RULE_SOURCE) != RULE_SHA:
        raise ValueError('registered existing state parser changed')
    spec = importlib.util.spec_from_file_location('E067_reviewed_rules', RULE_SOURCE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def single_spare_slot(body, *, terminal_flip=False):
    active = set()
    allowed = {DONE,LIVE,TARGET,OTHER_HELD,*DOWNSTREAM}
    for line in body.splitlines():
        parts = [v.strip() for v in line.split('|')]
        if len(parts) != 5:
            raise ValueError('malformed queue')
        job,name,state,partition,cpus = parts
        if not name.startswith('hfo2'):
            continue
        if job not in allowed or partition != 'hfacnormal01' or cpus != '32':
            raise ValueError('unexpected study allocation')
        if state in ('RUNNING','COMPLETING','CONFIGURING'):
            active.add(job)
        elif state != 'PENDING':
            raise ValueError('unexpected study state')
    expected = set() if terminal_flip else {LIVE}
    if active != expected:
        raise ValueError('registered live/terminal flip and at most one used slot required')
    return sorted(active)


def verify_terminal_receipt(receipt):
    result = receipt.get('terminal_observation', {})
    force = result.get('replayed_fmax_eV_A')
    if (receipt.get('job_id') != DONE or receipt.get('ordinary_converged') is not True
            or receipt.get('status') != 'audited_ordinary_converged_not_TS_or_sampling_certificate'
            or receipt.get('new_DFT_calls') != 0 or receipt.get('scheduler_mutations') != 0
            or receipt.get('physical_inputs_or_running_source_changed') is not False
            or receipt.get('holdout_generated_or_read') is not False
            or receipt.get('all_calls_original_six_physical_bytes_native_DSIZE32_and_raw_EFS_checked') is not True
            or result.get('ordinary_residual_pass') is not True
            or isinstance(force, bool) or not isinstance(force, (int,float))
            or not math.isfinite(force) or not 0 <= force <= .10):
        raise ValueError('actual converged terminal evidence required before a rolling handoff')


def pending_seed_preflight(seed):
    import numpy as np
    from ase.io import read
    import examples.hfo2_fixed_input_factory as transport
    from scripts.audit_hfo2_static_replica import audited_results
    from vcneb import validate_path_geometry, validate_periodic_path_lift
    manifest = json.loads((seed/'manifest.json').read_text())
    params = json.loads((seed/'factory_parameters.json').read_text())
    if manifest['physical_contract_sha256'] != transport.CONTRACT:
        raise ValueError('registered pending physical contract differs')
    source = Path(params['source_directory'])
    if any(sha(source/n) != h for n,h in transport.CONTRACT.items()):
        raise ValueError('pending original six physical inputs changed')
    images = read(seed/'seed.traj', index=':')
    if len(images) != 9 or any(a.calc is not None for a in images):
        raise ValueError('nine fresh calculator-free starter geometries required')
    validate_periodic_path_lift(images)
    validate_path_geometry(images, cell_scale=manifest['cell_scale_A'], minimum_distance=1.6, maximum_deformation=.25)
    if any(a.get_chemical_symbols() != ['Hf']*4+['O']*8 or not np.allclose(
            a.cell.array[:2], np.asarray(manifest['reference_cell_A'])[:2], atol=1e-10, rtol=0) for a in images):
        raise ValueError('pending ordered atoms or common substrate changed')
    caches = params['seed_cache_records']
    if len(caches) != 9 or any(r is not None for r in caches[1:-1]):
        raise ValueError('seven fresh moving images and two cached endpoints required')
    endpoint_audits = []
    for index in (0,8):
        cached = caches[index]
        raw = Path(cached['directory'])
        if (any(sha(raw/n) != h for n,h in cached['input_sha256'].items())
                or any(cached['input_sha256'][n] != h for n,h in transport.CONTRACT.items())
                or sha(raw/'OUT.ABACUS/running_scf.log') != cached['raw_log_sha256']
                or not transport.same_ordered_geometry(images[index], transport.read_fixed_hfo2_stru(raw/'STRU'))):
            raise ValueError('pending exact ordered endpoint raw cache changed')
        result = audited_results(raw)
        if (not np.isfinite(result['energy']) or result['forces'].shape != (12,3)
                or result['stress'].shape != (6,)):
            raise ValueError('pending full native endpoint EFS incomplete')
        endpoint_audits.append(dict(image_index=index, raw_source=str(raw),
            original_six_physical_bytes_native_DSIZE32_full_EFS_checked=True,
            raw_log_sha256=cached['raw_log_sha256']))
    return dict(geometry_images_checked=9, fresh_interiors=7, cached_endpoints=endpoint_audits, new_DFT_calls=0)


def release_once(audit_sha256):
    journal_path = ROOT/'rolling_handoff_journal.jsonl'
    if journal_path.exists():
        raise FileExistsError('inspect previous intent; never repeat a write or submission')
    audit_path = ROOT/'audit_receipt.json'
    if sha(audit_path) != audit_sha256:
        raise ValueError('reviewed terminal audit changed')
    receipt = json.loads(audit_path.read_text())
    verify_terminal_receipt(receipt)
    work = Path(receipt['workdir'])
    if sha(work/'vcneb_summary.json') != receipt['terminal_summary_sha256']:
        raise ValueError('terminal summary changed')
    rules = load_rules()
    parent = rules.run(['sacct','-X','-n','-P','-j',DONE,
        '--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,User'])
    if parent != DONE+'|COMPLETED|0:0|32|hfacnormal01|iai806':
        raise ValueError('actual parent completion not corroborated')
    flip = rules.run(['sacct','-X','-n','-P','-j',LIVE,
        '--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,User'])
    if flip not in {LIVE+'|'+s+'|0:0|32|hfacnormal01|iai806' for s in ('RUNNING','COMPLETING','COMPLETED')}:
        raise ValueError('registered flip accounting not corroborated; do not infer from missing squeue')
    terminal_flip = flip.split('|')[1] == 'COMPLETED'
    preflight_path = WAIT_ROOT/'preflight.json'
    if sha(preflight_path) != PREFLIGHT_SHA:
        raise ValueError('original pending-task preflight changed')
    preflight = json.loads(preflight_path.read_text())
    source = Path(preflight['source_root'])
    if (sha(source/rules.PILOT.rsplit('/source/', 1)[1]) != preflight['pilot_sha256']
            or sha(rules.BINARY) != preflight['ABACUS_binary_sha256']
            or any(sha(source/n) != h for n,h in preflight['source_code_sha256'].items())):
        raise ValueError('pending immutable runtime or binary changed')
    seed_record = next(r for r in preflight['reports'] if r['channel'] == 'PO_to_T')
    seed = Path(seed_record['seed_root'])
    if (sha(seed/'manifest.json') != seed_record['manifest_sha256']
            or any(sha(seed/n) != h for n,h in seed_record['seed_files_sha256'].items())
            or (WAIT_ROOT/'PO_to_T/band').exists()):
        raise ValueError('existing pending seed changed or production already exists')
    pending = pending_seed_preflight(seed)
    required = set() if terminal_flip else {LIVE}
    original = {j: rules.held_waiter(rules.run(['scontrol','show','job','-o',j]), j,
        required=required, allowed={DONE,LIVE}) for j in (TARGET,OTHER_HELD)}
    active = single_spare_slot(rules.run(['squeue','-h','-u','iai806','-o','%i|%100j|%T|%P|%C']), terminal_flip=terminal_flip)
    with journal_path.open('x') as journal:
        def record(event, **data):
            journal.write(json.dumps(dict(event=event, check_CST=datetime.now().isoformat(), **data))+'\n')
            journal.flush()
            os.fsync(journal.fileno())
        dependency = 'afterok:'+DONE
        record('rolling_handoff_verified', job_id=TARGET, flip_accounting=flip,
               terminal_flip=terminal_flip, audit_receipt_sha256=audit_sha256,
               active_jobs=active, held_rows=original, pending_seed_preflight=pending)
        try:
            if not terminal_flip:
                record('dependency_update_intent', job_id=TARGET, dependency=dependency)
                rules.run(['scontrol','update','JobId='+TARGET,'Dependency='+dependency])
                deadline = time.monotonic()+45
                while True:
                    body = rules.run(['scontrol','show','job','-o',TARGET])
                    try:
                        verified = rules.held_waiter(body, TARGET, required=set(), allowed={DONE})
                        break
                    except ValueError:
                        if time.monotonic() >= deadline:
                            raise TimeoutError('accepted dependency write unconfirmed; reconcile, never repeat')
                        time.sleep(4)
                record('dependency_independently_read_verified', job_id=TARGET, row=verified)
            single_spare_slot(rules.run(['squeue','-h','-u','iai806','-o','%i|%100j|%T|%P|%C']), terminal_flip=terminal_flip)
            rules.held_waiter(rules.run(['scontrol','show','job','-o',OTHER_HELD]), OTHER_HELD,
                              required=required, allowed={DONE,LIVE})
            record('existing_target_release_intent', job_id=TARGET)
            rules.run(['scontrol','release',TARGET])
            record('existing_target_released', job_id=TARGET)
        except Exception as exc:
            record('reconcile_existing_writes_do_not_repeat', exception=type(exc).__name__, message=str(exc))
            raise
    result = dict(status='existing_zero_strain_PO_to_T_released_into_single_spare_slot',
        check_CST=datetime.now().isoformat(), target_job=TARGET, audited_completed_parent=DONE,
        registered_flip_job=LIVE, flip_accounting_state=flip.split('|')[1],
        other_successor_still_held=OTHER_HELD,
        new_submissions=0, dependency_writes=0 if terminal_flip else 1, release_writes=1,
        max_study_simultaneous_allocations=2, new_independent_chains=0,
        physical_inputs_or_source_changed=False, holdout_generated_or_read=False,
        execution_started_claimed=False, executed_helper_sha256=sha(Path(__file__)),
        terminal_audit_sha256=audit_sha256, journal_sha256=sha(journal_path))
    with (ROOT/'rolling_handoff_receipt.json').open('x') as stream:
        json.dump(result, stream, indent=2)
    print(json.dumps(result))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reviewed-audit-sha256', required=True)
    release_once(parser.parse_args().reviewed_audit_sha256)
