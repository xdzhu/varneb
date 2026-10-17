"""One read-only snapshot of the two original E268 paths; no restart or DFT."""
import argparse
from datetime import datetime
import json
from pathlib import Path
import subprocess

import numpy as np

from examples.hfo2_fixed_input_factory import CONTRACT
from scripts.audit_hfo2_static_replica import audited_results, sha256

ROOT = Path('/public/home/iai806/abacus/agent-runs/20261018-varneb-physical-E268-r1')
JOBS = ('28938772', '28938773')
CHANNELS = ('PO_flip_T_pattern_preserving', 'PO_to_T')


def inspect(prefix=False):
    handles = [json.loads((ROOT/f'path_accepted_{i}.json').read_bytes()) for i in (0, 1)]
    if tuple(h['job_id'] for h in handles) != JOBS:
        raise ValueError('only the two original accepted handles may be inspected')
    raw = subprocess.check_output(['sacct', '-n', '-X', '-P', '-j', ','.join(JOBS),
        '--format=JobIDRaw,State,ExitCode,Partition,AllocCPUS,AllocNodes,ElapsedRaw,Start,End,NodeList'], text=True)
    scheduler = {r.split('|')[0]: r.split('|') for r in raw.strip().splitlines()
                 if r.split('|')[0] in JOBS}
    rows = []
    for job, channel in zip(JOBS, CHANNELS):
        work = ROOT/'paths'/channel/'band'
        calls = sorted(work.glob('image_*/scf_*/call_audit.json'))
        attempted = sorted(work.glob('image_*/scf_*'))
        row = dict(job_id=job, channel=channel, completed_SCFs=len(calls), attempted_SCFs=len(attempted))
        log = work/'vcneb.opt.log'
        row['optimizer_tail'] = log.read_text().splitlines()[-3:] if log.exists() else []
        summary = work/'vcneb_summary.json'
        row['summary'] = json.loads(summary.read_bytes()) if summary.exists() else None
        failure = work/'vcneb_failure.json'
        row['failure'] = json.loads(failure.read_bytes()) if failure.exists() else None
        times = [json.loads(p.read_bytes())['elapsed_seconds'] for p in calls]
        row['completed_call_seconds'] = dict(min=min(times), median=float(np.median(times)), max=max(times)) if times else None
        output = ROOT/f'path-{job}.out'
        row['slurm_tail'] = output.read_text(errors='replace').splitlines()[-8:] if output.exists() else []
        if prefix and calls:
            p = calls[0]; call = json.loads(p.read_bytes()); folder = p.parent
            if any(call['input_sha256'][n] != h or sha256(folder/n) != h for n, h in CONTRACT.items()):
                raise ValueError('first completed interior changed the fixed physical inputs')
            if sha256(folder/'STRU') != call['input_sha256']['STRU'] or sha256(folder/'OUT.ABACUS/running_scf.log') != call['raw_log_sha256']:
                raise ValueError('native geometry or log differs from call evidence')
            results = audited_results(folder)
            if any(not np.allclose(results[k], call['results'][k], atol=1e-12, rtol=0) for k in results):
                raise ValueError('first native E/F/stress does not replay')
            row['first_completed_native_audit'] = dict(call=str(folder), audit_sha256=sha256(p),
                physical_sha256=CONTRACT, raw_log_sha256=call['raw_log_sha256'],
                electronic_converged=True, real_MPI_ranks=32, energy_force_stress_complete=True,
                elapsed_seconds=call['elapsed_seconds'], prefix_only=True)
        rows.append(row)
    return dict(observed_at=datetime.now().astimezone().isoformat(), scheduler=scheduler,
                channels=rows, new_DFT_calls=0, restarts=0)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--audit-prefix', action='store_true')
    p.add_argument('--output', type=Path)
    a = p.parse_args(); report = inspect(a.audit_prefix)
    if a.output:
        with a.output.open('x', encoding='utf-8') as f:
            json.dump(report, f, indent=2, allow_nan=False); f.write('\n')
    print(json.dumps(report, allow_nan=False))
