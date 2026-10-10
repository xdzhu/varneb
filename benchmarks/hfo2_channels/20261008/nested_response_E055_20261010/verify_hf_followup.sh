#!/bin/bash
# The first script stopped after an analytic pass: ICU has no pytest.
# Keep its receipt/source; no package installation or rerun/overwrite.
set -euo pipefail
analysis_root=/public/home/iai806/abacus/agent-runs/20261010-varneb-nested-controls-E055-r1
[[ -f "$analysis_root/analytic_hf.json" && ! -e "$analysis_root/prior_pilot_scf_costs.json" ]] || exit 2
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONPATH="$analysis_root/source"
cd "$analysis_root/source"
py=/public/home/iai806/.conda/envs/icu/bin/python
"$py" - "$analysis_root" <<'PY'
import importlib.util,json,sys
from pathlib import Path
import numpy as np
from vcneb import ResponseEvaluationCost,recorded_response_dataset_cost
root=Path(sys.argv[1])
analytic=json.loads((root/'analytic_hf.json').read_text())
assert analytic['paired_model_points']==48 and max(analytic['maximum_errors'].values())<1e-9
assert analytic['DFT_calls']==0 and not analytic['HfO2_advantage_proven']
assert importlib.util.find_spec('pytest') is None
known=ResponseEvaluationCost('synthetic_shared','a'*64,'completed',120.,32)
failed=ResponseEvaluationCost('synthetic_failed','b'*64,'failed',60.,32)
unknown=ResponseEvaluationCost('synthetic_unknown','c'*64,'completed',None,None)
p=recorded_response_dataset_cost([known,known,failed],required_evaluation_ids=['synthetic_shared'])
assert p['unique_source_DFT_evaluations']==2 and p['failed_evaluations']==1
assert abs(p['recorded_transport_cpu_core_hours_sum']-1.6)<1e-12
q=recorded_response_dataset_cost([known,unknown])
assert not q['cost_complete'] and q['recorded_transport_cpu_core_hours_sum'] is None
report=dict(status='analytic_and_cost_semantics_pass_no_HF_pytest',python=sys.version,numpy=np.__version__,
  analytic_points=48,pytest_available=False,local_clean_full_regression_recorded_separately=True,
  package_installed_or_upgraded=False,new_DFT_calls=0,holdout_generated_or_read=False,
  original_harness_failure='No module named pytest after independent analytic pass; source and outputs preserved')
with (root/'hf_followup_verification.json').open('x') as f:json.dump(report,f,indent=2);f.write('\n')
print(json.dumps(report))
PY
"$py" -m scripts.export_hfo2_scf_costs \
 --run-root /public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-G2-E053-r1 \
 --output "$analysis_root/prior_pilot_scf_costs.json"
# No live code, physical parameters, allocation or held-out data changed.
