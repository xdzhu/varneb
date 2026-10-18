"""Package only two geometric seeds and finite controls; never run DFT."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys
import tarfile

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[4]))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    with path.open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(value, handle, indent=2); handle.write('\n')


def package(case, output):
    from examples.hfo2_fixed_input_factory import CONTRACT
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True); bundle = output/'bundle'; bundle.mkdir()
    for name in ('control.py', 'paths.slurm', 'protocol.md'):
        # Windows checkout normalization must not change the reviewed script.
        data = (case/name).read_bytes()
        if b'\r' in data:
            raise ValueError('package scripts and protocol must use LF')
        shutil.copyfile(case/name, bundle/name)
    shutil.copytree(case/'seeds', bundle/'seeds', ignore=shutil.ignore_patterns('__pycache__'))
    save(bundle/'config.json', dict(channels=['PO_flip_T_pattern_preserving', 'PO_to_T'],
        source_strain=.01, target_strain=.005, physical_input_sha256=CONTRACT,
        runtime_root='/public/home/iai806/abacus/agent-runs/20261018-varneb-physical-E268-r1',
        runtime_bundle_sha256='689a26fed0f47dede6d7c2f2230d3181fe5a16bbc5eee0356da6f033d89987cf'))
    files = {p.relative_to(bundle).as_posix(): sha(p) for p in sorted(bundle.rglob('*')) if p.is_file()}
    save(bundle/'bundle_files.json', files)
    archive = output/'bundle.tar.gz'
    with tarfile.open(archive, 'x:gz') as tar:
        tar.add(bundle, arcname='bundle')
    result = dict(bundle_sha256=sha(bundle/'bundle_files.json'), archive_sha256=sha(archive),
                  archive_bytes=archive.stat().st_size, files=len(files), new_DFT_calls=0)
    save(output/'package_receipt.json', result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(package(Path(__file__).resolve().parent, args.output)))
