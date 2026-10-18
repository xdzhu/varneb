"""Package only the reviewed two-band seed and controls; never submit."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import tarfile

CASE = Path(__file__).resolve().parent


def package(output):
    output = Path(output); output.mkdir(parents=True, exist_ok=False)
    bundle = output/'bundle'; bundle.mkdir()
    for name in ('control.py', 'paths.slurm', 'protocol.md'):
        if b'\r' in (CASE/name).read_bytes(): raise ValueError('runtime must use LF')
        shutil.copyfile(CASE/name, bundle/name)
    shutil.copytree(CASE/'seeds', bundle/'seeds')
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    files = {p.relative_to(bundle).as_posix(): sha(p) for p in sorted(bundle.rglob('*')) if p.is_file()}
    (bundle/'bundle_files.json').write_text(json.dumps(files, indent=2)+'\n', encoding='utf-8')
    archive = output/'bundle.tar.gz'
    with tarfile.open(archive, 'x:gz') as tar: tar.add(bundle, arcname='bundle')
    result = dict(bundle_sha256=sha(bundle/'bundle_files.json'), archive_sha256=sha(archive),
                  archive_bytes=archive.stat().st_size, files=len(files), new_DFT_calls=0)
    (output/'package_receipt.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    print(json.dumps(package(p.parse_args().output)))
