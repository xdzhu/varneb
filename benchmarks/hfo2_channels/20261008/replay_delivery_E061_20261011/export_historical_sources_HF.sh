#!/bin/bash
# Archive eight already executed E059 source files; no DFT or input reads.
set -euo pipefail
/public/home/iai806/.conda/envs/icu/bin/python - <<'PY'
from datetime import datetime
import hashlib
import json
from pathlib import Path
import zipfile

runs = Path('/public/home/iai806/abacus/agent-runs')
old = runs/'20261010-varneb-endpoint-work-E059-r1'
root = runs/'20261011-varneb-replay-delivery-E061-r1'
def sha(raw):return hashlib.sha256(raw).hexdigest()
receipt_raw = (old/'replay_receipt.json').read_bytes()
assert sha(receipt_raw) == 'cfa39ded1bf6fd09bc9ab77652ef85ecdc642322a81531f56dc4032633b6f431'
receipt = json.loads(receipt_raw)
expected = {'vcneb/strain_work.py','scripts/analyze_hfo2_endpoint_strain_work.py',
    'scripts/prepare_hfo2_clamped_chains.py','scripts/audit_hfo2_static_replica.py',
    'examples/hfo2_fixed_input_factory.py','vcneb/abacus.py',
    'paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.md','paper/VARNEB_JCTC/METHODS_DRAFT.md'}
assert set(receipt['source_sha256']) == expected
sources = {name:(old/'source'/name).read_bytes() for name in expected}
assert all(sha(raw) == receipt['source_sha256'][name] for name,raw in sources.items())
assert not root.exists(), 'single-use evidence export; reconcile existing receipt'
root.mkdir()
archive = root/'source_bundle.zip'
with archive.open('xb') as stream:
    with zipfile.ZipFile(stream,'w',compression=zipfile.ZIP_DEFLATED) as out:
        for name in sorted(sources):out.writestr(name,sources[name])
proof = dict(check_CST=datetime.now().isoformat(),historical_E059_receipt_sha256=sha(receipt_raw),
    historical_tested_tree=receipt['tested_staged_tree'],historical_source_archive_sha256=receipt['archive_sha256'],
    exact_raw_source_members=receipt['source_sha256'],source_bundle_sha256=sha(archive.read_bytes()),
    raw_source_bytes=sum(len(raw) for raw in sources.values()),source_files=8,
    new_numerical_HF_replay=False,new_DFT_calls=0,DFT_inputs_or_holdout_read=False,
    physical_parameters_or_running_sources_or_jobs_changed=False,
    environment_installation_or_HF_pytest=False,licensed_assets_exported=False,
    exact_historical_evidence_not_current_editable_manuscript_contract=True)
with (root/'source_export_receipt.json').open('x') as f:
    json.dump(proof,f,indent=2);f.write('\n')
print(json.dumps(dict(check_CST=proof['check_CST'],source_files=8,
    source_bundle_sha256=proof['source_bundle_sha256'],new_DFT_calls=0,new_numerical_HF_replay=False)))
PY
# Protect the heredoc terminator from Windows stdin's trailing CRLF.
