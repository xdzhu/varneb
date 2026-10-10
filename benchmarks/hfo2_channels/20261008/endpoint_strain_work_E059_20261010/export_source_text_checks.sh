#!/bin/bash
# Preserve exact HF archive identities and separately audit legacy CRLF text.
set -euo pipefail
/public/home/iai806/.conda/envs/icu/bin/python - <<'PY'
import hashlib,json
from pathlib import Path
root=Path('/public/home/iai806/abacus/agent-runs/20261010-varneb-endpoint-work-E059-r1')
receipt=json.loads((root/'replay_receipt.json').read_text());checks={}
for name,expected in receipt['source_sha256'].items():
 raw=(root/'source'/name).read_bytes()
 actual=hashlib.sha256(raw).hexdigest();assert actual==expected
 checks[name]=dict(raw_executed_sha256=actual,
  LF_canonical_text_sha256=hashlib.sha256(raw.replace(b'\r\n',b'\n')).hexdigest(),
  raw_CRLF_count=raw.count(b'\r\n'),raw_LF_count=raw.count(b'\n'))
report=dict(actual_source_files=8,raw_source_still_matches_original_replay=True,
 checks=checks,text_canonicalization_is_NOT_physical_input_normalization=True,
 scientific_sources_and_manuscript_require_exact_raw_bytes=True,
 source_files_modified=False,new_DFT_calls=0)
with (root/'source_text_checks.json').open('x') as f:json.dump(report,f,indent=2);f.write('\n')
print(json.dumps(dict(source_files=8,source_files_modified=False,new_DFT_calls=0)))
PY
# Canonical text is a comparison only; no source or physical file is rewritten.
