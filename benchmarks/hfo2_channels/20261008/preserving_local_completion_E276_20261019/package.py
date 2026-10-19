"""Package only committed runtime bytes, without networking or scheduler actions."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tarfile

CASE=Path(__file__).resolve().parent
REPO=CASE.parents[3]


def package(output):
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip()
    selected=[CASE/'worker.py',CASE/'run.slurm']+sorted(p for p in (CASE/'prepared').rglob('*') if p.is_file())
    payload={}
    for path in selected:
        name=path.relative_to(CASE).as_posix();relative=path.relative_to(REPO).as_posix()
        raw=subprocess.check_output(['git','show',commit+':'+relative],cwd=REPO)
        if raw!=path.read_bytes():raise ValueError('runtime file must exactly match committed bytes: '+name)
        payload[name]=raw
    files={name:hashlib.sha256(raw).hexdigest() for name,raw in payload.items()}
    manifest=(json.dumps(files,indent=2,sort_keys=True)+'\n').encode()
    payload['bundle_files.json']=manifest
    archive=output/'runtime.tar.gz'
    with tarfile.open(archive,'w:gz') as tar:
        for name,raw in payload.items():
            item=tarfile.TarInfo('bundle/'+name);item.size=len(raw);item.mtime=0;item.mode=0o644
            tar.addfile(item,io.BytesIO(raw))
    report=dict(commit=commit,files=files,bundle_pin=hashlib.sha256(manifest).hexdigest(),
        archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),archive_bytes=archive.stat().st_size,
        remote_root='/public/home/iai806/abacus/agent-runs/20261019-varneb-preserving-local-E276-r1',
        runtime_files=len(files),new_DFT_calls=0,scheduler_actions=0)
    with (output/'package.json').open('x',encoding='utf-8',newline='\n') as f:
        json.dump(report,f,indent=2);f.write('\n')
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    r=package(p.parse_args().output);print(json.dumps({k:v for k,v in r.items() if k!='files'}))
