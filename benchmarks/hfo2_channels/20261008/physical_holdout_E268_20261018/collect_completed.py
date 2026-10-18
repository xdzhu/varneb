"""Run the existing terminal readers on HF and export completed observations.

Execute through the frozen runtime's Python with its source on PYTHONPATH.
No calculator, optimizer or submission is invoked. Existing exports are refused.
"""
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import tarfile

ROOT = Path('/public/home/iai806/abacus/agent-runs/20261018-varneb-physical-E268-r1')
OUT = ROOT/'completed_paths_20261018'
ARCHIVE = ROOT/'completed_paths_20261018.tar.gz'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    if OUT.exists() or ARCHIVE.exists():
        raise FileExistsError('completed export already exists; inspect it before reuse')
    spec = importlib.util.spec_from_file_location('frozen_terminal_reader', Path(__file__).with_name('audit_paths.py'))
    reader = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reader)
    OUT.mkdir()
    results = []
    for index, short in enumerate(('preserving', 'T')):
        report = reader.audit(index, OUT/short)
        work = ROOT/'paths'/report['channel']/'band'
        shutil.copyfile(work/'vcneb.opt.log', OUT/short/'vcneb.opt.log')
        for i, call in enumerate(report['calls']):
            source = Path(call['directory'])
            destination = OUT/short/'all_calls'/f'{i:03d}'
            (destination/'OUT.ABACUS').mkdir(parents=True)
            for name in ('INPUT', 'KPT', 'STRU', 'call_audit.json', 'OUT.ABACUS/running_scf.log'):
                shutil.copyfile(source/name, destination/name)
        results.append({k: v for k, v in report.items() if k != 'calls'})
    files = {p.relative_to(OUT).as_posix(): sha(p) for p in sorted(OUT.rglob('*')) if p.is_file()}
    (OUT/'files.json').write_text(json.dumps(files, indent=2)+'\n', encoding='utf-8')
    with tarfile.open(ARCHIVE, 'x:gz') as archive:
        archive.add(OUT, arcname=OUT.name)
    print(json.dumps(dict(archive=str(ARCHIVE), archive_sha256=sha(ARCHIVE),
                         bytes=ARCHIVE.stat().st_size, files=len(files), results=results)))


if __name__ == '__main__':
    main()
