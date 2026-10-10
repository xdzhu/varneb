set -euo pipefail
runtime=/public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-G2-E054-r1
source_root="$runtime/source-fixed"
[[ "${FULL_REGRESSION_PASSED:-0}" == 1 && ! -e "$runtime/submission_handles.txt" ]] || exit 2
[[ "$(sha256sum "$runtime/source.tar" | cut -d' ' -f1)" == 48c7c0e8550e3812d53d810c5e9d720f03693b0a7bc085856a663afc2e7712c2 ]] || exit 2
[[ "$(sha256sum "$runtime/source-fixed.tar" | cut -d' ' -f1)" == 170170ed6c3f8749c4930ca8af6a67fc17689dcdf57481e916577e40aea70d9c ]] || exit 2
[[ "$(sha256sum "$source_root/cluster/hf_hfo2_clamped_chain_resume_20261010.slurm" | cut -d' ' -f1)" == 79d9990b42a8da90fe915550332d8060761d0ec108551035e284b69f7768ef02 ]] || exit 2
[[ "$(sha256sum "$runtime/full-regression.junit.xml" | cut -d' ' -f1)" == "${FULL_REGRESSION_JUNIT_SHA256:?set verified JUnit hash}" ]] || exit 2
for job in 28574708 28574709; do
 row=$(sacct -n -X -j "$job" --format=JobIDRaw,State,ExitCode -P | tr -d ' ')
 [[ "$row" == "$job|COMPLETED|0:0" ]] || exit 2
done
export PYTHONPATH="$source_root" PYTHONDONTWRITEBYTECODE=1
/public/home/iai806/.conda/envs/icu/bin/python - <<'PY'
import json
from pathlib import Path
import xml.etree.ElementTree as ET
from scripts.audit_hfo2_static_replica import sha256
runtime=Path('/public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-G2-E054-r1')
suites=ET.parse(runtime/'full-regression.junit.xml').getroot().findall('testsuite')
assert len(suites)==1 and suites[0].attrib['failures']==suites[0].attrib['errors']=='0'
assert int(suites[0].attrib['tests'])==1403 and int(suites[0].attrib['skipped'])==2
audit=json.loads((runtime/'audit_prepare.json').read_text())
assert audit['status']=='audited_clamped_geometry_continuations_prepared' and audit['new_DFT_calls']==0
assert not audit['holdout_generated_or_read'] and not audit['physical_inputs_changed']
verification=json.loads((runtime/'verification_after_transport_error.json').read_text())
assert verification['status']=='passed' and verification['new_DFT_calls']==0
assert verification['audit_prepare_sha256']==sha256(runtime/'audit_prepare.json')
for name in ('examples/hfo2_fixed_input_factory.py','scripts/export_hfo2_clamped_observation.py',
             'scripts/prepare_hfo2_clamped_resume.py','vcneb/core.py','vcneb/epitaxial_boundary.py',
             'cluster/hf_hfo2_clamped_chain_resume_20261010.slurm'):
    assert sha256(runtime/'source'/name)==sha256(runtime/'source-fixed'/name)
for r in audit['reports']:
    assert r['all_original_six_physical_bytes_and_native_DSIZE32_checked']
    assert r['current_frame_nine_exact_caches_checked'] and r['external_launch_prohibited']
    assert r['same_independent_chain'] and not r['FIRE_state_restored']
    seed=runtime/'seeds'/r['channel']
    assert sha256(seed/'manifest.json')==r['geometry_continuation_manifest_sha256']
    m=json.loads((seed/'manifest.json').read_text())
    assert all(sha256(seed/n)==h for n,h in m['files_sha256'].items())
PY
study_count=$(squeue -h -u iai806 -o '%j' | awk '$1 ~ /^hfo2/ {n++} END {print n+0}')
[[ "$study_count" == 0 ]] || { echo 'study allocation still present; do not duplicate' >&2; exit 2; }
for channel in PO_flip_T_pattern_preserving PO_to_M; do
 [[ ! -e "$runtime/$channel" ]] || exit 2
 job=$(sbatch --parsable --job-name="hfo2-G2-cont-${channel}" \
  --export="ALL,RUN_DFT=1,SOURCE_ROOT=$source_root,SEED_ROOT=$runtime/seeds/$channel,WORKDIR=$runtime/$channel/band" \
  --output="$runtime/${channel}.slurm.out" --error="$runtime/${channel}.slurm.err" \
  "$source_root/cluster/hf_hfo2_clamped_chain_resume_20261010.slurm")
 [[ "$job" =~ ^[0-9]+$ ]] || { echo 'ambiguous handle; reconcile scheduler, never repeat' >&2; exit 2; }
 printf '%s|%s\n' "$channel" "$job" | tee -a "$runtime/submission_handles.txt"
done
date +%FT%T%z
squeue -h -u iai806 -o '%i|%j|%T|%M|%C|%R'
# One reviewed continuation pair only; no recurring monitor or automatic retry.
