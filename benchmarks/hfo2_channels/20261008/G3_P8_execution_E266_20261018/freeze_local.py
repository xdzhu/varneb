"""Copy exact pre-label artifacts and freeze an inert local deployment; no DFT."""
from pathlib import Path
import hashlib
import importlib.util
import json

HERE=Path(__file__).resolve().parent;BASE=HERE.parent
spec=importlib.util.spec_from_file_location('E266_local_worker',HERE/'worker.py')
worker=importlib.util.module_from_spec(spec);spec.loader.exec_module(worker)
sources={
    'forecasts.json':BASE/'G3_anharmonic_model_E265_20261018/P8_forecasts.json',
    'starts.json':BASE/'G3_P_branch_independent_E263_20261018/starts_E264.json',
    'protocol.md':BASE/'G3_anharmonic_model_E265_20261018/P8_PROTOCOL.md'}
if any((HERE/name).exists() for name in (*sources,'config.json','bundle_files.json')):
    raise FileExistsError('fresh local artifact names required')
for name,path in sources.items():
    raw=path.read_bytes();key=name.split('.')[0]
    if hashlib.sha256(raw).hexdigest()!=worker.PINS[key]:raise ValueError('published pre-label artifact differs')
for name,path in sources.items():
    with (HERE/name).open('xb') as handle:handle.write(path.read_bytes())
old,_=worker.old_worker();old.study.exclusive_json(HERE/'config.json',worker.prepare(HERE))
files={name:worker.sha(HERE/name) for name in ('worker.py','run.slurm','config.json','forecasts.json','starts.json','protocol.md')}
old.study.exclusive_json(HERE/'bundle_files.json',files)
print(json.dumps(dict(config_sha256=worker.sha(HERE/'config.json'),bundle_manifest_sha256=worker.sha(HERE/'bundle_files.json'),new_DFT=0)))
