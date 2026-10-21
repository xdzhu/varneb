"""Package a fresh immutable sidecar; do not copy or modify the DFT runtime."""
import argparse
import json
from pathlib import Path
import shutil
import sys
import tarfile
REPO=Path(__file__).resolve().parents[4];sys.path.insert(0,str(REPO))
from scripts.audit_hfo2_static_replica import sha256


def package(output):
    case=Path(__file__).resolve().parent;output=Path(output)
    output.mkdir(parents=True,exist_ok=False);bundle=output/'bundle';bundle.mkdir()
    for name in ('worker.py','policy.py','run.slurm','protocol.md'):
        if b'\r' in (case/name).read_bytes():raise ValueError('execution files must use LF')
        shutil.copyfile(case/name,bundle/name)
    shutil.copyfile(case.with_name('preserving_local_completion_E276_20261019')/'worker.py',bundle/'native_start.py')
    shutil.copytree(case/'prepared',bundle/'prepared')
    files={p.relative_to(bundle).as_posix():sha256(p) for p in sorted(bundle.rglob('*')) if p.is_file()}
    with (bundle/'bundle_files.json').open('x',encoding='utf-8',newline='\n') as f:
        json.dump(files,f,indent=2);f.write('\n')
    with tarfile.open(output/'bundle.tar.gz','x:gz') as tar:tar.add(bundle,arcname='bundle')
    result=dict(bundle_pin=sha256(bundle/'bundle_files.json'),archive_sha256=sha256(output/'bundle.tar.gz'),files=len(files))
    with (output/'receipt.json').open('x',encoding='utf-8') as f:json.dump(result,f,indent=2)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    print(json.dumps(package(p.parse_args().output)))
