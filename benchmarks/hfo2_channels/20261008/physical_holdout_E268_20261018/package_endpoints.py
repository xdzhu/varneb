"""Create an inert source/bundle package; no SSH, scheduler or DFT calls."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import tarfile


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, obj):
    with path.open('x', encoding='utf-8', newline='\n') as f:
        json.dump(obj, f, indent=2); f.write('\n')


def package(case, source_archive, output):
    if output.exists(): raise FileExistsError(output)
    source_hashes = {}
    with tarfile.open(source_archive) as tar:
        for item in tar.getmembers():
            if item.isfile():
                name = Path(item.name)
                if name.is_absolute() or '..' in name.parts: raise ValueError('unsafe own archive path')
                source_hashes[item.name] = hashlib.sha256(tar.extractfile(item).read()).hexdigest()
    output.mkdir()
    bundle = output/'bundle'; bundle.mkdir()
    shutil.copyfile(source_archive, output/'source.tar')
    for name in ('endpoint_worker.py', 'endpoints.slurm', 'dispatch_endpoints.py', 'registration.json', 'protocol.md'):
        shutil.copyfile(case/name, bundle/name)
    shutil.copytree(case/'endpoint_seeds', bundle/'endpoint_seeds')
    save(bundle/'runtime.json', dict(source_archive_sha256=sha(source_archive), source_files_sha256=source_hashes))
    files = {p.relative_to(bundle).as_posix(): sha(p) for p in sorted(bundle.rglob('*')) if p.is_file()}
    save(bundle/'bundle_files.json', files)
    save(output/'package_receipt.json', dict(bundle_sha256=sha(bundle/'bundle_files.json'),
         source_archive_sha256=sha(source_archive), source_file_count=len(source_hashes),
         new_DFT_calls=0, target_interiors_generated=False))
    print((output/'package_receipt.json').read_text())


def path_package(case, forecasts, forecast_sha256, seeds, output):
    if output.exists(): raise FileExistsError(output)
    if sha(forecasts/'manifest.json') != forecast_sha256:
        raise ValueError('externally frozen forecast manifest differs')
    seed_report = json.loads((seeds/'chains.json').read_bytes())
    if seed_report['forecast_manifest_sha256'] != forecast_sha256 or len(seed_report['seeds']) != 2:
        raise ValueError('two forecast-bound seed chains required')
    output.mkdir(); bundle = output/'path_bundle'; bundle.mkdir()
    for name in ('path_control.py', 'paths.slurm'): shutil.copyfile(case/name, bundle/name)
    shutil.copytree(forecasts, bundle/'forecasts'); shutil.copytree(seeds, bundle/'seeds')
    save(bundle/'config.json', dict(forecast_manifest_sha256=forecast_sha256,
        channels=['PO_flip_T_pattern_preserving', 'PO_to_T']))
    files = {p.relative_to(bundle).as_posix(): sha(p) for p in sorted(bundle.rglob('*')) if p.is_file()}
    save(bundle/'bundle_files.json', files)
    save(output/'path_package_receipt.json', dict(path_bundle_sha256=sha(bundle/'bundle_files.json'),
         forecast_manifest_sha256=forecast_sha256, new_DFT_calls=0))
    print((output/'path_package_receipt.json').read_text())


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode', choices=('endpoints', 'paths'), default='endpoints')
    p.add_argument('--case', type=Path, default=Path(__file__).resolve().parent)
    p.add_argument('--source-archive', type=Path)
    p.add_argument('--forecasts', type=Path); p.add_argument('--forecast-sha256')
    p.add_argument('--seeds', type=Path)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if a.mode == 'endpoints': package(a.case, a.source_archive, a.output)
    else: path_package(a.case, a.forecasts, a.forecast_sha256, a.seeds, a.output)
