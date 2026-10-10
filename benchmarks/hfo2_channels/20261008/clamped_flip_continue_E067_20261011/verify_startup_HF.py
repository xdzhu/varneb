"""Independent read-only startup snapshot; never submit or modify a job."""
from datetime import datetime
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path('/public/home/iai806/abacus/agent-runs/20261011-varneb-clamped-flip-E067-r1')
JOBS = '28722810,28722320,28692775,28692776,28709788,28709789,28709790,28709791'


def run(argv):
    return subprocess.run(argv, check=True, capture_output=True, text=True, timeout=45).stdout.strip()


if __name__ == '__main__':
    work = ROOT/'PO_flip_T_pattern_preserving/band'
    log = work/'vcneb.opt.log'
    native = work/'image_0001/scf_000000/OUT.ABACUS/running_scf.log'
    report = dict(check_CST=datetime.now().isoformat(),
        raw_registered_queue=run(['squeue','-h','-j',JOBS,'-o','%i|%100j|%T|%P|%C']),
        raw_accounting=run(['sacct','-X','-n','-P','-j','28722810,28722320',
                            '--format=JobIDRaw,State,ExitCode,Elapsed,AllocCPUS,Start,Partition']),
        raw_waiters={j:run(['scontrol','show','job','-o',j]) for j in ('28692775','28692776')},
        new_flip_optimizer_log=log.read_text() if log.exists() else None,
        first_fresh_SCF_native_DSIZE_lines=[s for s in native.read_text().splitlines() if 'DSIZE' in s] if native.exists() else [],
        stderr_tail=(ROOT/'flip.slurm.err').read_text()[-2000:] if (ROOT/'flip.slurm.err').exists() else None,
        first_fresh_SCF_completion_claimed=False, new_scheduler_mutations=0,
        executed_verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    with (ROOT/'startup_verification.json').open('x') as stream:
        json.dump(report, stream, indent=2)
        stream.write('\n')
    print(json.dumps(report))
