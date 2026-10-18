"""Only completed force-evaluated frames from the original handle."""
import importlib.util
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/'benchmarks/hfo2_channels/20261008/preserving_connection_E273_20261019/collect_snapshot.py'
spec=importlib.util.spec_from_file_location('E273_snapshot_test',SCRIPT)
reader=importlib.util.module_from_spec(spec);spec.loader.exec_module(reader)


@pytest.mark.parametrize('state',['RUNNING','COMPLETED'])
def test_only_original_live_or_completed_resource_identity(state):
    raw=f'28949435|{state}|0:0|32|hfacnormal01|iai806|1|600|start|end'
    assert reader.require_original_handle(raw)[1]==state
    for wrong in (raw.replace('|32|','|64|'),raw.replace('0:0','1:0'),raw.replace('28949435','28948425')):
        with pytest.raises(ValueError):reader.require_original_handle(wrong)


@pytest.mark.parametrize('state',['PENDING','TIMEOUT','FAILED','CANCELLED'])
def test_other_states_are_not_promoted(state):
    with pytest.raises(ValueError):reader.require_original_handle(f'28949435|{state}|0:0|32|hfacnormal01|iai806|1|600|start|end')


def test_only_recorded_consecutive_bounded_evaluations_are_complete():
    text=' Step Time Energy fmax\nCheckedFIRE: 0 00:43:06 -9782.631831 0.125611\n'
    assert reader.completed_steps(text)==[0]
    assert reader.completed_steps(text+'CheckedFIRE: 1 00:54:01 -9782.632000 0.124000\n')==[0,1]
    for bad in ('',text+'CheckedFIRE: 2 00:55:01 -9782.632000 0.124000\n',
                text+'CheckedFIRE: 0 00:55:01 -9782.632000 0.124000\n'):
        with pytest.raises(ValueError):reader.completed_steps(bad)
