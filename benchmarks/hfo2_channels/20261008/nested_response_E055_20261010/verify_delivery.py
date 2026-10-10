"""Pin the actual E055 code, receipts, costs and retained harness failures."""
import argparse
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

from vcneb import ResponseEvaluationCost,recorded_response_dataset_cost


def run(case):
    case=Path(case)
    output=case/'validation.json'
    if output.exists():raise FileExistsError('fresh final validation only')
    def load(name):return json.loads((case/name).read_text())
    def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
    for name,field in (('analytic_local.json','source_sha256'),('analytic_hf.json','source_sha256'),
                       ('prior_pilot_scf_costs.json','export_source_sha256')):
        for source,sha in load(name)[field].items():
            assert digest(case/'executed_source'/source)==sha
    suite=ET.parse(case/'local_clean_full_junit.xml').getroot().find('testsuite')
    assert int(suite.attrib['tests'])==1461 and int(suite.attrib['failures'])==int(suite.attrib['errors'])==0
    cost=load('prior_pilot_scf_costs.json')
    records=[ResponseEvaluationCost(**r) for r in cost['records']]
    replay=recorded_response_dataset_cost(records,required_evaluation_ids=[r.evaluation_id for r in records])
    assert replay==cost['cost'] and replay['unique_source_DFT_evaluations']==154
    material=load('material_gate.json')
    assert not material['material_predictions_emitted'] and not material['held_out_structure_generated_or_labels_read']
    result=dict(status='verified_analytic_controls_and_prior_cost_not_material_forecast',
      tested_source_commit='86d2637ebcd6c502df43b91ca03977fd50120716',
      tested_source_archive_sha256='cbe6552a3fa4237fec2098c08c8056735a9a95df4e18691d09fed6b4f24a613c',
      clean_regression=dict(passed=1459,skipped=2,failed=0,errors=0,seconds=318.96,
          junit_sha256=digest(case/'local_clean_full_junit.xml')),
      focused_related_passed=175,local_delivery_passed=6,
      HF_numeric_points=48,HF_pytest_available=False,HF_package_installed_or_upgraded=False,
      completed_prior_SCFs=154,prior_SCF_transport_core_hours=replay['recorded_transport_cpu_core_hours_sum'],
      prior_scheduler_allocation_core_hours=cost['scheduler_allocation_cpu_core_hours'],
      material_gate_rejected_incomplete_observations=True,material_prediction_accuracy_measured=False,
      B2_to_B5_material_forecast_frozen=False,HfO2_advantage_proven=False,
      new_DFT_calls=0,new_independent_chains=0,physical_inputs_changed=False,
      held_out_structure_generated_or_labels_read=False,live_source_modified=False,
      new_CI_release_PyPI_or_automation=False,
      source_immutability='Exact executed-source subset preserved; no historical hash substituted with current checkout bytes',
      retained_zero_DFT_harness_issues=[
        'HF analytic checker passed, then pytest absent; separate no-install follow-up passed',
        'Direct local helper invocation omitted repository PYTHONPATH; explicit import path passed',
        'git diff used non-Git clean archive cwd; repeated in repository without mutation',
        'Initial progress snapshot checked optional wrong failure filename; corrected new receipt kept both',
        'First delivery replay5passed1failed on LF-working-vs-CRLF-executed parser; archived exact executed source,6passed'],
      files_sha256={p.relative_to(case).as_posix():digest(p) for p in sorted(case.rglob('*')) if p.is_file()},
      goal_complete=False,limitations=['Analytic algebra and cost replay are not prospective material accuracy',
        'Cost covers two completed first-pilot interiors ONLY, not total study/predictor preparation',
        'Live optimizer logs are not complete raw snapshot or TS audits; segment ETA is not convergence ETA'])
    with output.open('x',encoding='utf-8',newline='\n') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:result[k] for k in ('status','clean_regression','prior_SCF_transport_core_hours','goal_complete')}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',type=Path,required=True)
    args=parser.parse_args();run(args.case)
