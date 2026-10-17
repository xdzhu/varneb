"""One-shot held endpoint submission; ambiguity is reconciled, never retried."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess

ROOT = Path('/public/home/iai806/abacus/agent-runs/20261018-varneb-physical-E268-r1')


def command(args):
    return subprocess.check_output(args, text=True).strip()


def save(path, obj):
    with path.open('x', encoding='utf-8') as f:
        json.dump(obj, f, indent=2); f.write('\n')


def queue():
    raw = command(['squeue', '-h', '-u', 'iai806', '-o', '%i|%j|%T|%P|%C'])
    return [line for line in raw.splitlines() if line.split('|')[1].lower().startswith(('hfo2', 'varneb', 'vcneb'))]


def validate_held(raw, shard, job):
    fields = dict(re.findall(r'(\S+)=(\S+)', raw))
    expected = dict(JobId=job, JobState='PENDING', Reason='JobHeldUser', Partition='hfacnormal01',
                    NumCPUs='32', NumTasks='1', TimeLimit='08:00:00', JobName=f'hfo2-E268-endpoints-{shard}')
    expected['CPUs/Task'] = '32'
    if fields.get('NumNodes') not in ('1', '1-1') or any(fields.get(k) != v for k, v in expected.items()):
        raise ValueError('held resources differ; leave held')


def main(mode, pin):
    bundle = ROOT/'bundle'
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    if sha(bundle/'bundle_files.json') != pin:
        raise ValueError('bundle pin changed')
    files = json.loads((bundle/'bundle_files.json').read_bytes())
    if any(sha(bundle/n) != h for n, h in files.items()):
        raise ValueError('bundle drift')
    proof = json.loads((ROOT/'native_preflight.json').read_bytes())
    if proof['bundle_sha256'] != pin or len(proof['decks']) != 3:
        raise ValueError('three native preflight decks required')
    if mode == 'submit-held':
        if queue() or any((ROOT/f'{name}_{s}.json').exists() for s in (0, 1)
                          for name in ('intent', 'held', 'accepted', 'released')):
            raise ValueError('existing study allocation or submission intent; do not duplicate')
        for shard in (0, 1):
            if len(queue()) >= 2: raise ValueError('two allocation cap')
            cmd = ['sbatch', '--hold', '--parsable', f'--job-name=hfo2-E268-endpoints-{shard}',
                   f'--export=ALL,RUN_DFT=1,SHARD={shard},BUNDLE_SHA256={pin}',
                   f'--output={ROOT}/endpoints-%j.out', str(bundle/'endpoints.slurm')]
            save(ROOT/f'intent_{shard}.json', dict(command=cmd, bundle_sha256=pin))
            reply = command(cmd); job = reply.split(';')[0]
            if not job.isdigit(): raise ValueError('unknown scheduler reply; reconcile, never repeat')
            row = dict(job_id=job, shard=shard, bundle_sha256=pin, reply=reply,
                       raw_scontrol=command(['scontrol', 'show', 'job', '-o', job]))
            save(ROOT/f'held_{shard}.json', row); print(json.dumps(row), flush=True)
    else:
        rows = []
        if len(queue()) > 2: raise ValueError('two allocation cap')
        for shard in (0, 1):
            if (ROOT/f'accepted_{shard}.json').exists() or (ROOT/f'released_{shard}.json').exists():
                raise ValueError('original handle already accepted/released; inspect')
            row = json.loads((ROOT/f'held_{shard}.json').read_bytes())
            if row['bundle_sha256'] != pin: raise ValueError('held pin mismatch')
            raw = command(['scontrol', 'show', 'job', '-o', row['job_id']])
            validate_held(raw, shard, row['job_id'])
            rows.append((shard, row, raw))
        for shard, row, raw in rows:
            save(ROOT/f'accepted_{shard}.json', dict(**row, validated_scontrol=raw))
        for shard, row, _ in rows:
            save(ROOT/f'release_intent_{shard}.json', dict(job_id=row['job_id']))
            reply = command(['scontrol', 'release', row['job_id']])
            record = dict(job_id=row['job_id'], reply=reply,
                          raw_scontrol=command(['scontrol', 'show', 'job', '-o', row['job_id']]))
            save(ROOT/f'released_{shard}.json', record); print(json.dumps(record), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=('submit-held', 'release'))
    p.add_argument('--bundle-sha256', required=True)
    a = p.parse_args(); main(a.mode, a.bundle_sha256)
