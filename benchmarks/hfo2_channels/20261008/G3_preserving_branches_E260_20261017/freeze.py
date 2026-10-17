"""Assemble and freeze the small external bundle; no DFT or scheduler calls."""
from pathlib import Path
import importlib.util
import hashlib
import json
import shutil

CASE=Path(__file__).resolve().parent
ROOT=CASE.parents[3]
SOURCES={
    'legacy_soft_worker.py':CASE.parent/'G3_soft_probe_E256_20261016/worker.py',
    'branch_relaxation.py':ROOT/'vcneb/branch_relaxation.py',
    'soft_direction.json':CASE.parent/'G3_released_complete_E254_20261016/soft_direction.json',
    'soft_probe_analysis.json':CASE.parent/'G3_soft_probe_E256_20261016/analysis_E257.json'}


def freeze():
    for name,source in SOURCES.items():
        target=CASE/name
        if target.exists():
            if target.read_bytes()!=source.read_bytes():raise ValueError('copied source drift')
        else:shutil.copyfile(source,target)
    spec=importlib.util.spec_from_file_location('E260_freeze_worker',CASE/'worker.py')
    worker=importlib.util.module_from_spec(spec);spec.loader.exec_module(worker)
    worker.study.exclusive_json(CASE/'config.json',worker.prepare(CASE))
    names=['worker.py','branch_relaxation.py','legacy_soft_worker.py','run.slurm',
           'config.json','protocol.md','soft_direction.json','soft_probe_analysis.json']
    worker.study.exclusive_json(CASE/'bundle_files.json',
        {name:hashlib.sha256((CASE/name).read_bytes()).hexdigest() for name in names})
    print(json.dumps(dict(config_sha256=worker.sha(CASE/'config.json'),
        bundle_manifest_sha256=worker.sha(CASE/'bundle_files.json'),new_DFT=0)))


if __name__=='__main__':freeze()
