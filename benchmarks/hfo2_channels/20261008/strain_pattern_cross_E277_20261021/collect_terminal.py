"""Export the two original terminal jobs once; no new DFT or scheduler writes."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tarfile

ROOT=Path('/public/home/iai806/abacus/agent-runs/20261021-varneb-strain-pattern-E277-r1')
NAMES=('zero_q1','plus1_q0')
PIN='ecb498b40ec3f19b593bfca198a754d1f2cdc62deb389fcc1dbb6b0d85ea2683'


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def save(path,value):
    with Path(path).open('x',encoding='utf-8',newline='\n') as f:
        json.dump(value,f,indent=2,allow_nan=False);f.write('\n')


def collect(output):
    accepted={n:json.loads((ROOT/f'accepted_{n}.json').read_bytes()) for n in NAMES}
    rows=[]
    for name in NAMES:
        job=accepted[name]['job_id']
        raw=subprocess.check_output(['sacct','-n','-X','-P','-j',job,
            '--format=JobIDRaw,State,ExitCode,Partition,AllocNodes,AllocCPUS,ElapsedRaw,Start,End'],text=True)
        lines=[r.split('|')[:9] for r in raw.strip().splitlines()]
        if (len(lines)!=1 or lines[0][:6]!=[job,'COMPLETED','0:0','hfacnormal01','1','32']
                or accepted[name]['bundle_pin']!=PIN):
            raise ValueError(f'{name}: original job must have a verified normal terminal outcome: {raw}')
        complete=json.loads((ROOT/'results'/name/'complete.json').read_bytes())
        if complete['job_id']!=job or complete['bundle_pin']!=PIN:raise ValueError('terminal identity differs')
        rows.append(dict(condition=name,job_id=job,scheduler=lines[0],
            allocated_core_hours=float(lines[0][6])*32/3600,new_calls=complete['new_calls']))
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    native=output/'native_terminal';native.mkdir()
    for row in rows:
        name=row['condition'];source=ROOT/'results'/name;dest=native/name;dest.mkdir()
        shutil.copyfile(source/'complete.json',dest/'complete.json')
        records=list(source.glob('evaluations/eval-*/result.json'))
        if len(records)!=row['new_calls']:raise ValueError('all attempted evaluations must be accounted for')
        for record in records:
            result=json.loads(record.read_bytes())
            if result['status']!='complete':raise ValueError('failed electronic attempt requires explicit diagnosis')
            target=dest/record.parent.relative_to(source);(target/'scf_000000/OUT.ABACUS').mkdir(parents=True)
            for n in ('result.json','scf_000000/INPUT','scf_000000/KPT','scf_000000/STRU',
                      'scf_000000/call_audit.json','scf_000000/abacus.out','scf_000000/abacus.err',
                      'scf_000000/OUT.ABACUS/running_scf.log'):
                shutil.copyfile(record.parent/n,target/n)
    files={p.relative_to(native).as_posix():sha(p) for p in sorted(native.rglob('*')) if p.is_file()}
    save(native/'files.json',files)
    archive=output/'native_terminal.tar.gz'
    with tarfile.open(archive,'x:gz') as tar:tar.add(native,arcname='native_terminal')
    receipt=dict(rows=rows,bundle_pin=PIN,archive_sha256=sha(archive),files=len(files),
        collector_sha256=sha(__file__),new_DFT_calls=0,scheduler_mutations=0)
    save(output/'receipt.json',receipt);return receipt


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    print(json.dumps(collect(p.parse_args().output)))
