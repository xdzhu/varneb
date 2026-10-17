"""Synthetic admission records, not actual material completion evidence."""
import importlib.util
from pathlib import Path
import pytest

P=Path(__file__).resolve().parents[1]/'benchmarks/hfo2_channels/20261008/G3_preserving_branches_E260_20261017/audit.py'
s=importlib.util.spec_from_file_location('E260_audit_tests',P)
w=importlib.util.module_from_spec(s);s.loader.exec_module(w)
RAW='28934048|COMPLETED|0:0|3600|32|node3|hfacnormal01|1\n28934049|COMPLETED|0:0|1800|32|node28|hfacnormal01|1\n'


def test_allocation_cost_keeps_both_synthetic_shards():
    rows=w.accounting(RAW)
    assert [r['job_id'] for r in rows]==w.JOBS
    assert sum(r['allocation_core_hours'] for r in rows)==48


@pytest.mark.parametrize('raw',[
    RAW.replace('COMPLETED','RUNNING',1),RAW.replace('0:0','1:0',1),
    RAW.splitlines()[0],RAW.replace('|32|','|64|',1),RAW+RAW.splitlines()[0]+'\n'])
def test_nonterminal_partial_failure_and_wrong_resources_are_not_a_complete_audit(raw):
    with pytest.raises(ValueError):w.accounting(raw)
