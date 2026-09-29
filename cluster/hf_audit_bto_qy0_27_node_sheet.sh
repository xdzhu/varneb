#!/bin/bash
# Read-only raw-log audit of all 15 new Qy=0 conditional points.
# Run once, after both Slurm arrays have finished. Creates new audit JSON only.
set -euo pipefail
root=/public/home/iai806/abacus/agent-runs/20260929-varneb-bto-qy0-grid27
prior=/public/home/iai806/abacus/agent-runs/20260928-varneb-bto-qy0-pilot
code=${prior}/code-holdouts
inputs=${prior}/inputs
python=/public/home/iai806/.conda/envs/icu/bin/python
[[ -s "${root}/plan.json" ]] || exit 2
mapfile -t points < <("${python}" - "${root}/plan.json" <<'PY'
import json
from pathlib import Path
import sys
p = json.loads(Path(sys.argv[1]).read_text())
if p.get('status') != '15_missing_Qy0_inputs_preflighted_no_DFT' or len(p['new']) != 15:
    raise SystemExit('wrong conditional sheet plan')
for row in p['new']:
    print(row['name'])
PY
)
[[ "${#points[@]}" == 15 ]] || exit 2
cd "${code}"
for point in "${points[@]}"; do
  case "${point}" in
    qz01_qx02|qz02_qx01) result=${point}_stress_refined_result.json ;;
    *) result=${point}_result.json ;;
  esac
  work=${root}/run-${point}
  output=${root}/audit-${point}.json
  [[ -s "${work}/${result}" && ! -e "${output}" ]] || {
    echo "missing completed summary or audit already exists: ${point}" >&2; exit 3;
  }
  "${python}" scripts/audit_bto_q1q2_conditional_pilot.py \
    --workdir "${work}" --summary-name "${result}" \
    --preflight "${root}/preflights/${point}.json" \
    --report "${inputs}/bto_t_to_c_q1q2_strain_gradient_preflight_2026-09-25.json" \
    --reference "${inputs}/cubic_CONTCAR" \
    --force-constants "${inputs}/bto_cubic_gamma_force_constants.npz" \
    --phonopy-eigenpairs "${inputs}/bto_cubic_phonopy_gamma_eigenpairs.npz" \
    --gamma-provenance "${inputs}/bto_cubic_gamma_phonon_provenance.json" \
    --force-sets "${inputs}/bto_cubic_gamma_FORCE_SETS" \
    --eigenpairs-provenance "${inputs}/bto_cubic_phonopy_gamma_eigenpairs_provenance.json" \
    --grid-result-manifest "${inputs}/result_manifest_hf_27777454.json" \
    --output "${output}"
done
