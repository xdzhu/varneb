"""Historical HF source evidence plus current numerical replay, not byte lock-in."""
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

CASE = Path(__file__).resolve().parents[1]/'benchmarks/hfo2_channels/20261008/prediction_scoring_E064_20261011'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_actual_HF_receipt_pins_historical_numeric_files_not_current_editable_files():
    receipt = json.loads((CASE/'HF_receipt.json').read_text())
    for native, local in (('prediction_evaluation.py', 'executed_numeric.py'),
                           ('check_prediction_evaluation.py', 'executed_final_check.py')):
        assert sha(CASE/local) == receipt['source_sha256'][native]
    assert sha(CASE/'run_HF.py') == receipt['executed_wrapper_sha256']
    assert sha(CASE/'HF_numeric_check.json') == receipt['result_sha256']
    assert receipt['isolated_numeric_kernel_only']
    assert not receipt['full_HF_pytest_or_package_installation']
    assert not receipt['training_or_holdout_read']
    assert not receipt['running_source_jobs_inputs_changed']
    assert not receipt['material_forecast_frozen'] and not receipt['material_advantage_proven']


def test_current_numeric_implementation_replays_complete_historical_synthetic_panel():
    spec = importlib.util.spec_from_file_location('E064_historical_checker', CASE/'executed_final_check.py')
    checker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checker)  # historical inputs; current numeric implementation
    replay = checker.check()
    native = json.loads((CASE/'HF_numeric_check.json').read_text())
    entries = 0

    def compare(current, old, key=''):
        nonlocal entries
        if isinstance(old, dict):
            for child, value in old.items():
                if child == 'limitations':
                    continue  # editable explanations are not a numeric historical byte contract
                assert child in current
                compare(current[child], value, child)
        elif isinstance(old, list):
            assert len(current) == len(old)
            for a, b in zip(current, old): compare(a, b, key)
        elif isinstance(old, bool):
            assert current is old; entries += 1
        elif isinstance(old, (float, int)):
            assert current == pytest.approx(old, abs=1e-12, rel=0); entries += 1
        elif old is None:
            assert current is None
        elif key in ('model_id', 'case_id', 'energy_unit', 'comparison_contract_sha256',
                      'frozen_model_sha256', 'baseline_models', 'candidate_model'):
            assert current == old

    compare(replay, native)
    assert entries > 100
