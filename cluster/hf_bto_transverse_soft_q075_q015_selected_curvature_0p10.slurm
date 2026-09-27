#!/bin/bash
# Independent second-step fixed-Q Hessian for the selected +Q_y BTO holdout.
# Keep the five-atom Gamma source and ABACUS PBE/100 Ry/full 10 au DZP,
# electronic 4x4x4 mesh. This does not certify a conditional PES or barrier.
#SBATCH -J varneb-bto-q075-hess2
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=02:00:00
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260927-varneb-bto-soft-gridpilot-qx030/slurm-%x-%j.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260927-varneb-bto-soft-gridpilot-qx030/slurm-%x-%j.err

set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 ]] || exit 2
[[ "${SLURM_JOB_PARTITION:-}" == hfacnormal01 && "${SLURM_JOB_NUM_NODES:-0}" == 1 \
   && "${SLURM_NTASKS:-0}" == 1 && "${SLURM_CPUS_PER_TASK:-0}" -ge 32 ]] || exit 2

root=/public/home/iai806/abacus/agent-runs/20260927-varneb-bto-soft-gridpilot-qx030
code=${root}/code
inputs=${root}/inputs
workdir=${root}/run-q075_q015
parent=${workdir}/conditional_q075_q015_result.json
replay=${root}/audit-q075_q015-branch-replay-general-27787487.json
prior=${workdir}/conditional_q075_q015_selected_curvature_0p05.json
prior_raw_audit=${root}/audit-q075_q015-curvature-27787714.json
prior_curvature_audit=${root}/audit-q075_q015-curvature-reconstruction-27787714.json
geometry_preflight=${inputs}/bto_transverse_soft_q075_q015_selected_curvature_0p10_preflight_2026-09-27.json
result=${workdir}/conditional_q075_q015_selected_curvature_0p10.json
python=/public/home/iai806/.conda/envs/icu/bin/python
abacus=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus
pseudo_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Pseudopotential
basis_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Orb-DZP-10au
env_script=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh

[[ -s "${workdir}/contract.json" && -s "${parent}" && -s "${replay}" \
   && -s "${prior}" && -s "${prior_raw_audit}" \
   && -s "${prior_curvature_audit}" && -s "${geometry_preflight}" \
   && ! -e "${result}" ]] || exit 3
echo '82a3f36ca29db2118d9e0afbf88eb3cd76aaa5c5c15fd4547a0e3f0bdedd8f27  '"${parent}" | sha256sum -c -
echo '933bcd6fb9325e198ec2bc0236304c061bcd29de540c8e637cc63915386d1114  '"${replay}" | sha256sum -c -
echo '31874afd95b7b56b4a22ddf7f92cba712c2633682c0bda6b0e19c1aee3a26f08  '"${prior}" | sha256sum -c -
echo '62cb8d01091f4d72de50283318d8632a04c5f52a117a43e49bd40b762cf2a1de  '"${prior_raw_audit}" | sha256sum -c -
echo '925ad71e732b3c1faed2852ccbc128ca54d1f7fc1c63bb5d99e15cd9ae2a7a7c  '"${prior_curvature_audit}" | sha256sum -c -
echo '1589a26a00c5b5d775be8af267ae00775e26192966acff7ddbe84b6c635fc6f0  '"${geometry_preflight}" | sha256sum -c -
echo 'de0c1c3aeac376b9e4f4ea1f9ff742e85e70f9ba08823824c46c358d547a17bf  '"${code}/examples/run_bto_q1q2_conditional_abacus.py" | sha256sum -c -
echo '91b70d48a488b5a26f1e563bbce1ba0e75581a01d9145158d57108b8dac4bbfc  '"${code}/vcneb/mode_surface.py" | sha256sum -c -
"${python}" - "${parent}" "${prior}" "${prior_raw_audit}" "${prior_curvature_audit}" "${geometry_preflight}" <<'PY'
import hashlib, json, pathlib, sys
parent, prior, raw_audit, audit, preflight = map(pathlib.Path, sys.argv[1:])
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
s, c, r, a, p = (json.loads(path.read_text()) for path in
                 (parent, prior, raw_audit, audit, preflight))
if not (s['status'] == 'orthogonal_gradient_and_stress_converged_curvature_unchecked'
        and c['status'] == 'orthogonal_curvature_screen_complete_not_branch_certification'
        and c['start_from_result_sha256'] == sha(parent)
        and c['n_new_evaluations'] == 32
        and r['n_individually_audited_DFT_points'] == 109
        and r['source_sha256']['summary'] == sha(parent)
        and a['status'] == 'no_negative_curvature_at_one_difference_step_only'
        and a['source_sha256']['curvature_result'] == sha(prior)
        and a['source_sha256']['all_points_audit'] == sha(raw_audit)
        and a['n_matched_signed_probes'] == 32
        and p['kind'] == 'bto_transverse_soft_curvature_probe_geometry_preflight_no_dft'
        and p['status'] == 'all_signed_probes_geometry_safe_not_a_curvature_result'
        and p['source_sha256']['summary'] == sha(parent)
        and p['n_signed_probes'] == 32 and p['curvature_step_sqrt_amu_A'] == 0.10):
    raise SystemExit('second-step source or geometry chain changed')
PY

module purge
source "${env_script}"
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS I_MPI_HYDRA_BOOTSTRAP
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "${code}"
"${python}" examples/run_bto_q1q2_conditional_abacus.py \
  --preflight "${inputs}/bto_transverse_soft_conditional_q075_q015_preflight_2026-09-27.json" \
  --report "${inputs}/bto_t_to_c_q1q2_strain_gradient_preflight_2026-09-25.json" \
  --reference "${inputs}/cubic_CONTCAR" \
  --force-constants "${inputs}/bto_cubic_gamma_force_constants.npz" \
  --phonopy-eigenpairs "${inputs}/bto_cubic_phonopy_gamma_eigenpairs.npz" \
  --gamma-provenance "${inputs}/bto_cubic_gamma_phonon_provenance.json" \
  --force-sets "${inputs}/bto_cubic_gamma_FORCE_SETS" \
  --eigenpairs-provenance "${inputs}/bto_cubic_phonopy_gamma_eigenpairs_provenance.json" \
  --grid-result-manifest "${inputs}/result_manifest_hf_27777454.json" \
  --workdir "${workdir}" --pseudo-dir "${pseudo_dir}" --basis-dir "${basis_dir}" \
  --abacus-bin "${abacus}" --mpi-ranks 32 \
  --gradient-tolerance 0.003 --orthogonal-amplitude-bound 3.0 \
  --max-iterations 1 --max-new-evaluations 32 \
  --max-absolute-stress-kbar 1.0 --resume --run-dft \
  --start-from-result "${parent}" --curvature-only --curvature-step 0.10 \
  --result-name "$(basename "${result}")"
