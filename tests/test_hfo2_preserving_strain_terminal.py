"""Bounded terminal policy only; real jobs are never polled or launched here."""
import importlib.util
from pathlib import Path
import json
import pytest

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT/'benchmarks/hfo2_channels/20261008/preserving_strain_E275_20261019'


def collector():
    spec = importlib.util.spec_from_file_location('E275_collect_test', CASE/'collect_terminal.py')
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


def test_only_exact_original_completed_allocation_is_terminal():
    m = collector(); good = '28951226|COMPLETED|0:0|32|hfacnormal01|iai806|1|100|start|end'
    assert m.terminal_row(good, '28951226')[0] == '28951226'
    for raw in (good.replace('COMPLETED', 'RUNNING'), good.replace('COMPLETED', 'TIMEOUT'),
                good.replace('0:0', '1:0'), good.replace('|32|', '|1|'), good.replace('|1|100', '|2|100'), ''):
        with pytest.raises(ValueError): m.terminal_row(raw, '28951226')


def test_stepzero_uses_seven_fresh_calls_and_cap_is_not_convergence():
    m = collector(); row = m.terminal_row('28951226|COMPLETED|0:0|32|hfacnormal01|iai806|1|100|s|e', '28951226')
    params = json.loads((CASE/'seeds/strain_0050/factory_parameters.json').read_bytes())
    preflight = dict(factory=m.FACTORY, calculator_parameters=params)
    s = dict(final_max_generalized_force_eV_per_A=.09, climbing_image_active_final=False,
             status='converged', converged=True, termination='force_threshold')
    assert m.admit('28951226', row, s, 0, False, preflight, 'test-producer') == 7
    s.update(final_max_generalized_force_eV_per_A=.12, status='max_steps_reached',
             converged=False, termination='step_limit')
    assert m.admit('28951226', row, s, 60, False, preflight, 'test-producer') == 427
    for step in (59, 61):
        with pytest.raises(ValueError): m.admit('28951226', row, s, step, False, preflight, 'test-producer')
    with pytest.raises(ValueError): m.admit('28951226', row, s, 60, True, preflight, 'test-producer')


def test_fresh_reader_has_no_launch_or_cancellation_and_copies_preflight():
    text = (CASE/'collect_terminal.py').read_text()
    assert "'sacct'" in text and 'vcneb_preflight.json' in text
    for forbidden in ('sbatch', 'scancel', 'scontrol', 'material_runner', 'subprocess.run('):
        assert forbidden not in text
