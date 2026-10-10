"""Separate post-write snapshot; read Slurm only, never hold/release/submit."""
from datetime import datetime
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path('/public/home/iai806/abacus/agent-runs/20261011-varneb-terminal-guard-E066-r1')
JOBS = '28661019,28722320,28692775,28692776,28709788,28709789,28709790,28709791'


def run(argv):
    return subprocess.run(argv, check=True, capture_output=True, text=True, timeout=45).stdout.strip()


if __name__ == '__main__':
    rows = {job: run(['scontrol', 'show', 'job', '-o', job]) for job in ('28692775', '28692776')}
    report = dict(check_CST=datetime.now().isoformat(), raw_waiter_rows=rows,
        raw_registered_queue=run(['squeue', '-h', '-j', JOBS, '-o', '%i|%100j|%T|%P|%C']),
        raw_parent_accounting=run(['sacct', '-X', '-n', '-P', '-j', '28661019,28722320',
                                  '--format=JobIDRaw,State,ExitCode,Elapsed,AllocCPUS,Partition']),
        executed_verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        new_scheduler_mutations=0, new_DFT_calls=0)
    with (ROOT/'independent_verification.json').open('x') as stream:
        json.dump(report, stream, indent=2)
        stream.write('\n')
    print(json.dumps(report))
