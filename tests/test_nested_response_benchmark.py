"""Independent analytic benchmark, not a material-prediction score."""
import pytest

from scripts.check_nested_response import run_benchmark
from scripts.export_hfo2_scf_costs import parse_scheduler,export_costs


def test_bounded_independent_quadratic_checks():
    result=run_benchmark()
    assert result['paired_model_points']==48 and result['groups']==4
    assert max(result['maximum_errors'].values())<1e-9
    assert result['partial_measurement_statuses'][0]=='available_restricted_quadratic_only'
    assert all(s=='unavailable_unmeasured_curvature' for s in result['partial_measurement_statuses'][1:])
    assert result['training_instability_promoted_dimensions']==1
    assert result['DFT_calls']==0 and not result['HfO2_advantage_proven']


ROWS='28574708|COMPLETED|0:0|32|9564|hfacnormal01\n28574709|COMPLETED|0:0|32|7966|hfacnormal01\n'


def test_actual_scheduler_schema():
    parsed=parse_scheduler(ROWS)
    assert parsed['28574708']['allocated_cpu_cores']==32
    assert parsed['28574709']['allocation_seconds']==7966


@pytest.mark.parametrize('rows',[ROWS.replace('COMPLETED','RUNNING'),ROWS.replace('0:0','1:0'),
    ROWS.replace('|32|','|1|'),ROWS.replace('hfacnormal01','wrong_queue'),ROWS.splitlines()[0],
    ROWS+ROWS,ROWS.replace('9564','-1'),ROWS.replace('28574708','unknown')])
def test_scheduler_rejects_incomplete_failed_or_wrong_resource_proof(rows):
    with pytest.raises(ValueError):parse_scheduler(rows)


def test_cost_export_cannot_read_live_resume_or_other_namespace(tmp_path):
    with pytest.raises(ValueError,match='registered'):
        export_costs(tmp_path,ROWS)
