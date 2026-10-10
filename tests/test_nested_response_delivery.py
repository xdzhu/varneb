"""Replay the E055 evidence without remote/licensed assets or calculators."""
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import pytest

from vcneb import ResponseEvaluationCost,recorded_response_dataset_cost


ROOT=Path(__file__).resolve().parents[1]
CASE=ROOT/'benchmarks/hfo2_channels/20261008/nested_response_E055_20261010'


def load(name):
    return json.loads((CASE/name).read_text(encoding='utf-8'))


@pytest.mark.parametrize('name',['analytic_local.json','analytic_hf.json'])
def test_actual_numeric_receipts_are_analytic_not_material_accuracy(name):
    receipt=load(name)
    assert receipt['paired_model_points']==48 and max(receipt['maximum_errors'].values())<1e-9
    assert receipt['maximum_errors']['stable_B4_B5_equivalence']==0
    assert receipt['DFT_calls']==0 and receipt['calculator_calls']==0
    assert not receipt['HfO2_advantage_proven'] and not receipt['B2_to_B5_material_forecast_frozen']
    for p,h in receipt['source_sha256'].items():
        assert hashlib.sha256((CASE/'executed_source'/p).read_bytes()).hexdigest()==h


def test_actual154cost_records_reproduce_both_prior_receipts_and_no_fake_free_cache():
    receipt=load('prior_pilot_scf_costs.json')
    records=[ResponseEvaluationCost(**r) for r in receipt['records']]
    report=recorded_response_dataset_cost(records,required_evaluation_ids=[r.evaluation_id for r in records])
    assert report==receipt['cost']
    assert report['unique_source_DFT_evaluations']==154 and report['failed_evaluations']==0
    assert report['recorded_transport_cpu_core_hours_sum']==pytest.approx(154.45101636884112,abs=1e-10)
    assert receipt['scheduler_allocation_cpu_core_hours']==pytest.approx(155.82222222222222,abs=1e-10)
    prior=load('../clamped_G2_E054_20261010/audit_prepare.json')
    assert report['recorded_wall_seconds_sum']==pytest.approx(sum(r['prior_SCF_transport_seconds_sum'] for r in prior['reports']))
    reused=recorded_response_dataset_cost(records+records,required_evaluation_ids=[r.evaluation_id for r in records])
    assert reused==report  # duplicated cache consumers are not154newDFTcalls
    assert receipt['new_DFT_calls']==0 and not receipt['holdout_generated_or_read']
    assert len(receipt['per_call_proofs'])==154
    assert all(p['native_MPI_ranks']==32 and p['native_complete_EFS_checked']
               and p['six_original_physical_inputs_checked'] for p in receipt['per_call_proofs'])
    for p,h in receipt['export_source_sha256'].items():
        assert hashlib.sha256((CASE/'executed_source'/p).read_bytes()).hexdigest()==h


def test_failed_remote_harness_does_not_become_a_claimed_HF_pytest_pass():
    r=load('hf_followup_verification.json')
    assert not r['pytest_available'] and r['analytic_points']==48
    assert not r['package_installed_or_upgraded'] and r['new_DFT_calls']==0
    assert 'No module named pytest' in r['original_harness_failure']
    suite=ET.parse(CASE/'local_clean_full_junit.xml').getroot().find('testsuite')
    assert int(suite.attrib['tests'])==1461 and int(suite.attrib['skipped'])==2
    assert int(suite.attrib['failures'])==int(suite.attrib['errors'])==0


def test_actual_incomplete_material_observations_refused_before_prediction_output():
    r=load('material_gate.json')
    assert r['status']=='actual_pilot_observations_rejected_as_complete_training'
    assert r['ordinary_residual_pass']==[False,False] and all(f>.10 for f in r['fmax_eV_A'])
    assert not r['complete_G2_training_network_provided'] and not r['material_predictions_emitted']
    assert not r['held_out_structure_generated_or_labels_read'] and r['new_DFT_calls']==0
    for path,sha in r['input_sha256'].items():
        assert hashlib.sha256((ROOT/path.replace('\\','/')).read_bytes()).hexdigest()==sha
    assert not (CASE/'prediction_must_not_exist').exists()


def test_live_progress_observation_does_not_claim_raw_audit_or_mutate_jobs():
    r=load('live_jobs_snapshot_corrected.json')
    assert len(r['jobs'])==2 and r['new_DFT_calls']==0 and not r['jobs_mutated_or_submitted']
    assert not r['live_source_overwritten'] and not r['holdout_generated_or_read']
    assert r['production_script_sha256']=='79d9990b42a8da90fe915550332d8060761d0ec108551035e284b69f7768ef02'
    for job in r['jobs']:
        assert job['current_force_log_is_not_full_raw_snapshot_audit']
        assert job['estimate_is_segment_cap_not_convergence']
        assert job['failure_report_path'].endswith('/vcneb_failure.json')
        assert job['scheduler_row'][1]=='RUNNING' and not job['failure_report_present']


def test_validation_pins_exact_delivered_evidence_without_claiming_goal_complete():
    r=load('validation.json')
    for p,h in r['files_sha256'].items():
        assert hashlib.sha256((CASE/p).read_bytes()).hexdigest()==h
    assert r['new_DFT_calls']==0 and not r['physical_inputs_changed']
    assert not r['material_prediction_accuracy_measured'] and not r['goal_complete']
