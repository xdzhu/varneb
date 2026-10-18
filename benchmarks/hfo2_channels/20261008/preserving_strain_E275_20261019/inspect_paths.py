"""Read only the two accepted E275 handles; optionally reparse first native EFS."""
from datetime import datetime
import json
from pathlib import Path
import subprocess
import numpy as np
from examples.hfo2_fixed_input_factory import CONTRACT
from scripts.audit_hfo2_static_replica import audited_results, sha256

ROOT = Path('/public/home/iai806/abacus/agent-runs/20261019-varneb-preserving-strain-E275-r1')
PIN = '00ea6f9f82e50aeea60cfd560ecdcc76cec68b3308a0683e8bdf0fe922dfe8e3'


def inspect():
    rows = []
    for i, label in enumerate(('strain_0050', 'strain_0100')):
        accepted = json.loads((ROOT/f'accepted_{i}.json').read_bytes())
        if accepted['label'] != label or accepted['bundle_sha256'] != PIN:
            raise ValueError('only the original accepted two study handles may be inspected')
        job = accepted['job_id']; work = ROOT/'paths'/label/'band'
        raw = subprocess.check_output(['sacct', '-n', '-X', '-P', '-j', job,
            '--format=JobIDRaw,State,ExitCode,Partition,AllocCPUS,AllocNodes,ElapsedRaw,Start,End,NodeList'], text=True)
        calls = sorted(work.glob('image_*/scf_*/call_audit.json'))
        times = [json.loads(p.read_bytes())['elapsed_seconds'] for p in calls]
        row = dict(job_id=job, label=label, scheduler=raw.strip(), completed_SCFs=len(calls),
            optimizer_tail=(work/'vcneb.opt.log').read_text().splitlines()[-3:]
                if (work/'vcneb.opt.log').exists() else [],
            summary=json.loads((work/'vcneb_summary.json').read_bytes())
                if (work/'vcneb_summary.json').exists() else None,
            failure=json.loads((work/'vcneb_failure.json').read_bytes())
                if (work/'vcneb_failure.json').exists() else None,
            output_tail=(ROOT/f'path-{job}.out').read_text(errors='replace').splitlines()[-8:]
                if (ROOT/f'path-{job}.out').exists() else [],
            call_seconds_median=float(np.median(times)) if times else None)
        if calls:
            call = json.loads(calls[0].read_bytes()); folder = calls[0].parent
            if (any(call['input_sha256'][n] != digest or sha256(folder/n) != digest
                    for n, digest in CONTRACT.items())
                    or sha256(folder/'STRU') != call['input_sha256']['STRU']
                    or sha256(folder/'OUT.ABACUS/running_scf.log') != call['raw_log_sha256']):
                raise ValueError('first native input/geometry/log changed')
            parsed = audited_results(folder)
            if any(not np.allclose(parsed[k], call['results'][k], atol=1e-12, rtol=0) for k in parsed):
                raise ValueError('native first-call E/F/stress does not replay')
            row['first_native_EFS_replayed'] = dict(directory=str(folder), MPI_ranks=32,
                elapsed_seconds=call['elapsed_seconds'], input_sha256=call['input_sha256'],
                raw_log_sha256=call['raw_log_sha256'], prefix_not_full_path=True)
        rows.append(row)
    return dict(observed_at=datetime.now().astimezone().isoformat(), channels=rows,
                new_DFT_calls=0, scheduler_mutations=0, runtime_source_unchanged=True)


if __name__ == '__main__': print(json.dumps(inspect(), allow_nan=False))
