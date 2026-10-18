"""Future terminal admission uses native convergence, not exit-zero alone."""
import importlib.util
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
CASE=ROOT/'benchmarks/hfo2_channels/20261008/preserving_connection_E273_20261019'
spec=importlib.util.spec_from_file_location('E273_terminal_test',CASE/'collect_terminal.py')
reader=importlib.util.module_from_spec(spec);spec.loader.exec_module(reader)


def preflight():
    return dict(factory='examples.hfo2_fixed_input_factory:make_clamped_resume_cached_factory',
        calculator_parameters=dict(source_directory='unchanged',seed_cache_records=[{'pinned':'cache'} for _ in range(9)]))


def summary(converged):
    return dict(final_max_generalized_force_eV_per_A=.099 if converged else .12,
        climbing_image_active_final=False,status='converged' if converged else 'max_steps_reached',
        converged=converged,termination='force_threshold' if converged else 'step_limit')


def test_original_resource_accounting_and_no_live_or_timeout_promotion():
    raw='28949435|COMPLETED|0:0|32|hfacnormal01|iai806|1|600|start|end'
    assert reader.terminal_row(raw)[7]=='600'
    for state in ('RUNNING','PENDING','TIMEOUT','FAILED','CANCELLED'):
        with pytest.raises(ValueError):reader.terminal_row(raw.replace('COMPLETED',state))


def test_genuine_convergence_and_bounded_unconverged_outcome_remain_distinct():
    row=reader.terminal_row('28949435|COMPLETED|0:0|32|hfacnormal01|iai806|1|600|start|end')
    assert reader.admit(row,summary(True),5,False,preflight())==35
    assert reader.admit(row,summary(False),24,False,preflight())==168
    with pytest.raises(ValueError):reader.admit(row,summary(False),5,False,preflight())
    with pytest.raises(ValueError):reader.admit(row,summary(True),25,False,preflight())


def test_no_failure_or_unconverged_force_can_be_labelled_converged():
    row=reader.terminal_row('28949435|COMPLETED|0:0|32|hfacnormal01|iai806|1|600|start|end')
    with pytest.raises(ValueError):reader.admit(row,summary(True),5,True,preflight())
    bad=summary(True);bad['final_max_generalized_force_eV_per_A']=.101
    with pytest.raises(ValueError):reader.admit(row,bad,5,False,preflight())
    bad=summary(True);bad['climbing_image_active_final']=True
    with pytest.raises(ValueError):reader.admit(row,bad,5,False,preflight())


def test_initial_cache_policy_cannot_turn_into_seven_unrecorded_fresh_calls():
    row=reader.terminal_row('28949435|COMPLETED|0:0|32|hfacnormal01|iai806|1|600|start|end')
    broken=preflight();broken['calculator_parameters']['seed_cache_records'][4]=None
    with pytest.raises(ValueError):reader.admit(row,summary(True),5,False,broken)
