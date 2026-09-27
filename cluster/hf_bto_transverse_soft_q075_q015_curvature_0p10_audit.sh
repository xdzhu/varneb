#!/bin/bash
# Read-only second-step audit after Slurm 27787746 is terminal and successful.
set -euo pipefail

root=/public/home/iai806/abacus/agent-runs/20260927-varneb-bto-soft-gridpilot-qx030
code=${root}/code
inputs=${root}/inputs
workdir=${root}/run-q075_q015
analysis=${root}/analysis-code
python=/public/home/iai806/.conda/envs/icu/bin/python
parent=${workdir}/conditional_q075_q015_result.json
result=${workdir}/conditional_q075_q015_selected_curvature_0p10.json
canary=${root}/audit-q075_q015-canary-27787471.json
replay=${root}/audit-q075_q015-branch-replay-general-27787487.json
prior_raw_audit=${root}/audit-q075_q015-curvature-27787714.json
prior_curvature_audit=${root}/audit-q075_q015-curvature-reconstruction-27787714.json
raw_audit=${root}/audit-q075_q015-curvature-0p10-27787746.json
curvature_audit=${root}/audit-q075_q015-curvature-reconstruction-0p10-27787746.json

[[ -s "${result}" && -s "${parent}" && -s "${canary}" && -s "${replay}" \
   && -s "${prior_raw_audit}" && -s "${prior_curvature_audit}" \
   && ! -e "${raw_audit}" && ! -e "${curvature_audit}" ]] || exit 3
[[ "$(sacct -j 27787746 -X --format=State,ExitCode -n -P | tr -d ' ')" == 'COMPLETED|0:0' ]] || {
  echo 'Slurm job 27787746 is not successful and terminal' >&2; exit 4;
}
echo '82a3f36ca29db2118d9e0afbf88eb3cd76aaa5c5c15fd4547a0e3f0bdedd8f27  '"${parent}" | sha256sum -c -
echo '933bcd6fb9325e198ec2bc0236304c061bcd29de540c8e637cc63915386d1114  '"${replay}" | sha256sum -c -
echo '62cb8d01091f4d72de50283318d8632a04c5f52a117a43e49bd40b762cf2a1de  '"${prior_raw_audit}" | sha256sum -c -
echo '925ad71e732b3c1faed2852ccbc128ca54d1f7fc1c63bb5d99e15cd9ae2a7a7c  '"${prior_curvature_audit}" | sha256sum -c -
echo '8038e1d93e0a9e960e891ceff155aeddf39d79584269d46382807605ab28b3d2  '"${analysis}/audit_bto_q1q2_conditional_pilot_current.py" | sha256sum -c -
echo 'c6c49b38bbe5207e1564b30cabab2619ac5a7f19a8885f653919ac8a133204d8  '"${analysis}/audit_bto_q1q2_curvature_cache_current.py" | sha256sum -c -

export PYTHONPATH=${code}
"${python}" "${analysis}/audit_bto_q1q2_conditional_pilot_current.py" \
  --workdir "${workdir}" \
  --preflight "${inputs}/bto_transverse_soft_conditional_q075_q015_preflight_2026-09-27.json" \
  --report "${inputs}/bto_t_to_c_q1q2_strain_gradient_preflight_2026-09-25.json" \
  --reference "${inputs}/cubic_CONTCAR" \
  --force-constants "${inputs}/bto_cubic_gamma_force_constants.npz" \
  --phonopy-eigenpairs "${inputs}/bto_cubic_phonopy_gamma_eigenpairs.npz" \
  --gamma-provenance "${inputs}/bto_cubic_gamma_phonon_provenance.json" \
  --force-sets "${inputs}/bto_cubic_gamma_FORCE_SETS" \
  --eigenpairs-provenance "${inputs}/bto_cubic_phonopy_gamma_eigenpairs_provenance.json" \
  --canary-audit "${canary}" \
  --grid-result-manifest "${inputs}/result_manifest_hf_27777454.json" \
  --summary-name "$(basename "${parent}")" \
  --output "${raw_audit}"

"${python}" "${analysis}/audit_bto_q1q2_curvature_cache_current.py" \
  --workdir "${workdir}" \
  --preflight "${inputs}/bto_transverse_soft_conditional_q075_q015_preflight_2026-09-27.json" \
  --report "${inputs}/bto_t_to_c_q1q2_strain_gradient_preflight_2026-09-25.json" \
  --reference "${inputs}/cubic_CONTCAR" \
  --force-constants "${inputs}/bto_cubic_gamma_force_constants.npz" \
  --phonopy-eigenpairs "${inputs}/bto_cubic_phonopy_gamma_eigenpairs.npz" \
  --refined-result "${parent}" \
  --curvature-result "${result}" \
  --all-points-audit "${raw_audit}" \
  --branch-replay "${replay}" \
  --output "${curvature_audit}"
