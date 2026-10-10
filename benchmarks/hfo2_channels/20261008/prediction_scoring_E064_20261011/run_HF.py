"""Isolated NumPy-kernel check; not HF package installation/full pytest."""
from datetime import datetime
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import types

import numpy as np

ROOT = Path('/public/home/iai806/abacus/agent-runs/20261011-varneb-prediction-scoring-E064-r1')
PINNED = {
    'prediction_evaluation.py': 'b5bce40436b6956b59e049ab3fa5ccdd721fc84678be083a327c5b2a20a401f4',
    'check_prediction_evaluation.py': '3ecf00e0bdcc3a61b68c39ef4230d20ece531a944d60abd14508ca822de5864e',
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def main():
    if (ROOT/'HF_numeric_check.json').exists() or (ROOT/'HF_receipt.json').exists():
        raise FileExistsError('single-use isolated check; never overwrite evidence')
    if any(sha(ROOT/name) != value for name, value in PINNED.items()):
        raise ValueError('transported numeric/checker source differs')
    # Explicit isolation: no installed/production vcneb is imported and no
    # calculator, live directory, training structure or held-out label is read.
    namespace = types.ModuleType('vcneb')
    namespace.__path__ = []
    sys.modules['vcneb'] = namespace
    load('vcneb.prediction_evaluation', ROOT/'prediction_evaluation.py')
    checker = load('E064_check', ROOT/'check_prediction_evaluation.py')
    result = checker.check()
    with (ROOT/'HF_numeric_check.json').open('x') as stream:
        json.dump(result, stream, indent=2); stream.write('\n')
    receipt = dict(check_CST=datetime.now().isoformat(), status=result['status'],
        isolated_numeric_kernel_only=True, full_HF_pytest_or_package_installation=False,
        python=sys.version, numpy=np.__version__, source_sha256=PINNED,
        executed_wrapper_sha256=sha(Path(__file__)),
        result_sha256=sha(ROOT/'HF_numeric_check.json'), new_DFT_calls=0, new_job_submissions=0,
        running_source_jobs_inputs_changed=False, training_or_holdout_read=False,
        material_forecast_frozen=False, material_advantage_proven=False)
    with (ROOT/'HF_receipt.json').open('x') as stream:
        json.dump(receipt, stream, indent=2); stream.write('\n')
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
