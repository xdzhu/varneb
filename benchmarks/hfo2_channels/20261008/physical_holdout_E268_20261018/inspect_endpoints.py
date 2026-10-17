"""Read only exact original endpoint handles; optionally export small evidence."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile

ROOT = Path('/public/home/iai806/abacus/agent-runs/20261018-varneb-physical-E268-r1')
PHASES = ('PO_plus', 'PO_minus_T_preserving', 'T')


def inspect(complete=False):
    handles = [json.loads((ROOT/f'accepted_{s}.json').read_bytes()) for s in (0, 1)]
    jobs = [h['job_id'] for h in handles]
    if jobs != ['28938150', '28938151']:
        raise ValueError('only original endpoint handles')
    raw = subprocess.check_output(['sacct', '-n', '-X', '-j', ','.join(jobs), '-P',
        '--format=JobIDRaw,State,ExitCode,Partition,AllocCPUS,AllocNodes,ElapsedRaw,Start,End'], text=True)
    scheduler = {}
    for line in raw.splitlines():
        fields = line.strip().split('|')
        if fields[0] in jobs: scheduler[fields[0]] = fields
    states = []
    for phase in PHASES:
        folder = ROOT/'endpoints'/phase
        files = sorted((folder/'calculator/image_0000').glob('scf_*/call_audit.json'))
        state = folder/'endpoint_relax_state.json'
        summary = folder/'endpoint_relax_summary.json'
        row = dict(phase=phase, completed_SCFs=len(files),
                   attempted_SCFs=len(list((folder/'calculator/image_0000').glob('scf_*'))))
        path = summary if summary.exists() else state
        if path.exists():
            data = json.loads(path.read_bytes())
            row.update({k: data[k] for k in ('status', 'optimizer_steps',
                'max_atomic_force_eV_per_A', 'open_traction_norm_kbar')})
        if files:
            row['last_call_seconds'] = json.loads(files[-1].read_bytes())['elapsed_seconds']
        states.append(row)
    report = dict(scheduler=scheduler, phases=states, target_interiors_read=False, new_DFT_calls=0)
    if complete:
        if set(scheduler) != set(jobs) or any(row[1:6] != ['COMPLETED', '0:0', 'hfacnormal01', '32', '1']
                                             for row in scheduler.values()):
            raise ValueError('both original endpoint jobs must be completed0:0 on prescribed resources')
        paths = []
        for phase in PHASES:
            folder = ROOT/'endpoints'/phase
            audit = json.loads((folder/'native_audit.json').read_bytes())
            if not audit['endpoint_converged'] or audit['new_SCF_calls'] != next(r['completed_SCFs'] for r in states if r['phase'] == phase):
                raise ValueError('complete native converged endpoint audit required')
            paths += [folder/n for n in ('CONTCAR', 'endpoint_relax_summary.json', 'native_audit.json', 'relax.log')]
            for call in sorted((folder/'calculator/image_0000').glob('scf_*')):
                paths += [call/n for n in ('INPUT', 'KPT', 'STRU', 'call_audit.json',
                                          'OUT.ABACUS/running_scf.log', 'abacus.out', 'abacus.err')]
        report['allocation_core_hours'] = sum(float(row[6]) for row in scheduler.values())*32/3600
        report['exported_file_sha256'] = {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
        archive = ROOT/'endpoint_export_E268.tar'
        with tarfile.open(archive, 'x') as tar:
            for path in paths: tar.add(path, arcname=path.relative_to(ROOT).as_posix(), recursive=False)
        report['archive_sha256'] = hashlib.sha256(archive.read_bytes()).hexdigest()
        with (ROOT/'endpoint_export_E268.json').open('x', encoding='utf-8') as f:
            json.dump(report, f, indent=2); f.write('\n')
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--complete', action='store_true')
    a = p.parse_args(); print(json.dumps(inspect(a.complete)))
