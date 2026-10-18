"""One read-only snapshot of E269's accepted jobs, optionally native prefix.

The filename deliberately does not shadow Python's standard-library inspect.
"""
import argparse
from datetime import datetime
import json
from pathlib import Path
import subprocess
import sys

if __package__ in (None, ''):
    candidate = Path(__file__).resolve().parents[4]
    if (candidate/'examples/hfo2_fixed_input_factory.py').is_file():
        sys.path.insert(0, str(candidate))

import numpy as np

from examples.hfo2_fixed_input_factory import CONTRACT
from scripts.audit_hfo2_static_replica import audited_results, sha256

ROOT = Path('/public/home/iai806/abacus/agent-runs/20261018-varneb-midpoint-reverse-E269-r1')
CHANNELS = ('PO_flip_T_pattern_preserving', 'PO_to_T')


def inspect(prefix=False):
    accepted = [json.loads((ROOT/f'accepted_{i}.json').read_bytes()) for i in (0, 1)]
    jobs = [r['job_id'] for r in accepted]
    if len(set(jobs)) != 2 or [r['channel'] for r in accepted] != list(CHANNELS):
        raise ValueError('only the exact two accepted E269 jobs may be read')
    raw = subprocess.check_output(['sacct', '-n', '-X', '-P', '-j', ','.join(jobs),
        '--format=JobIDRaw,State,ExitCode,Partition,AllocCPUS,AllocNodes,ElapsedRaw,Start,End,NodeList'], text=True)
    scheduler = {r.split('|')[0]: r.split('|') for r in raw.strip().splitlines()
                 if r.split('|')[0] in jobs}
    rows = []
    for job, channel in zip(jobs, CHANNELS):
        work = ROOT/'paths'/channel/'band'
        calls = sorted(work.glob('image_*/scf_*/call_audit.json'))
        times = [json.loads(p.read_bytes())['elapsed_seconds'] for p in calls]
        log = work/'vcneb.opt.log'
        summary, failure = work/'vcneb_summary.json', work/'vcneb_failure.json'
        output = ROOT/f'path-{job}.out'
        row = dict(job_id=job, channel=channel, completed_SCFs=len(calls),
            attempted_SCFs=len(list(work.glob('image_*/scf_*'))),
            optimizer_tail=log.read_text().splitlines()[-3:] if log.exists() else [],
            summary=json.loads(summary.read_bytes()) if summary.exists() else None,
            failure=json.loads(failure.read_bytes()) if failure.exists() else None,
            slurm_tail=output.read_text(errors='replace').splitlines()[-8:] if output.exists() else [],
            completed_call_seconds=dict(min=min(times), median=float(np.median(times)), max=max(times)) if times else None)
        if prefix and calls:
            path = calls[0]; call = json.loads(path.read_bytes()); folder = path.parent
            if (any(call['input_sha256'][n] != h or sha256(folder/n) != h for n, h in CONTRACT.items())
                    or sha256(folder/'STRU') != call['input_sha256']['STRU']
                    or sha256(folder/'OUT.ABACUS/running_scf.log') != call['raw_log_sha256']):
                raise ValueError('native first-call physics, structure or output changed')
            parsed = audited_results(folder)
            if any(not np.allclose(parsed[k], call['results'][k], atol=1e-12, rtol=0) for k in parsed):
                raise ValueError('native first-call E/F/stress differs')
            row['first_completed_native_audit'] = dict(directory=str(folder),
                call_audit_sha256=sha256(path), raw_log_sha256=call['raw_log_sha256'],
                electronic_converged=True, genuine_MPI_ranks=32, complete_EFS=True,
                elapsed_seconds=call['elapsed_seconds'], prefix_only=True)
        rows.append(row)
    return dict(observed_at=datetime.now().astimezone().isoformat(), scheduler=scheduler,
                channels=rows, new_DFT_calls=0, scheduler_mutations=0)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit-prefix', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = inspect(args.audit_prefix)
    if args.output:
        with args.output.open('x', encoding='utf-8') as handle:
            json.dump(result, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps(result, allow_nan=False))
