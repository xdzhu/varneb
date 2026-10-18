"""Package the one-chain protocol and exact-cache seed; never submit or run."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import tarfile


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, data):
    with path.open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(data, handle, indent=2, allow_nan=False); handle.write('\n')


def package(case, output):
    output.mkdir(parents=True, exist_ok=False); bundle = output/'bundle'; bundle.mkdir()
    for name in ('control.py', 'run.slurm', 'protocol.md'):
        if b'\r' in (case/name).read_bytes():
            raise ValueError('runtime source and protocol must use LF')
        shutil.copyfile(case/name, bundle/name)
    shutil.copytree(case/'prepared', bundle/'prepared')
    files = {p.relative_to(bundle).as_posix(): sha(p) for p in sorted(bundle.rglob('*')) if p.is_file()}
    save(bundle/'bundle_files.json', files)
    archive = output/'bundle.tar.gz'
    with tarfile.open(archive, 'x:gz') as tar: tar.add(bundle, arcname='bundle')
    result = dict(bundle_pin=sha(bundle/'bundle_files.json'), archive_sha256=sha(archive),
                  files=len(files), archive_bytes=archive.stat().st_size, new_DFT_calls=0, scheduler_mutations=0)
    save(output/'package_receipt.json', result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); print(json.dumps(package(Path(__file__).resolve().parent, args.output)))
