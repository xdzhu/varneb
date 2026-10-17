"""Use an actual archived Slurm schema, without scheduler writes."""
import importlib.util
import json
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
CASE=ROOT/'benchmarks/hfo2_channels/20261008'
spec=importlib.util.spec_from_file_location('E260_dispatch_test',
    CASE/'G3_preserving_branches_E260_20261017/dispatch.py')
w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w)
RAW=json.loads((CASE/'G3_T_independent_E258_20261017/held_submission_0.json').read_bytes())['raw_scontrol']
RAW=RAW.replace('hfo2-G3-E258-T0','hfo2-G3-E260-P0')


@pytest.mark.parametrize('nodes',['1','1-1'])
def test_actual_slurm_cpu_per_task_spelling_and_exact_node_range(nodes):
    assert w.validate_held(RAW.replace('NumNodes=1-1','NumNodes='+nodes),0,'28913247')['CPUs/Task']=='32'


@pytest.mark.parametrize('old,new',[
    ('NumNodes=1-1','NumNodes=1-2'),('NumCPUs=32','NumCPUs=64'),
    ('CPUs/Task=32','CPUs/Task=16'),('Partition=hfacnormal01','Partition=other'),
    ('JobState=PENDING','JobState=RUNNING'),('TimeLimit=04:00:00','TimeLimit=24:00:00')])
def test_wrong_resource_or_state_is_refused(old,new):
    with pytest.raises(ValueError,match='held resources'):w.validate_held(RAW.replace(old,new),0,'28913247')
