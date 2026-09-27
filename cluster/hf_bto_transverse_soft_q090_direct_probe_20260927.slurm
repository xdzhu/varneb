#!/usr/bin/env bash
# Four signed lowest-direction statics per audited Qz=0.9 branch; no optimizer.
# ABACUS PBE/100 Ry/Ba-Ti-O 10 au DZP, electronic 4x4x4, 32 MPI.
# The Gamma force-constant source is the five-atom 1x1x1 cell.
#SBATCH -J varneb-bto-q090-soft4
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --array=0-1%2
#SBATCH --time=00:45:00
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260927-varneb-bto-soft-gridpilot-qx030/slurm-soft4-%A_%a.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260927-varneb-bto-soft-gridpilot-qx030/slurm-soft4-%A_%a.err

set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 ]] || { echo 'set RUN_DFT=1 explicitly' >&2; exit 2; }
[[ "${SLURM_JOB_PARTITION:-}" == hfacnormal01 && "${SLURM_JOB_NUM_NODES:-0}" == 1 \
   && "${SLURM_NTASKS:-0}" == 1 && "${SLURM_CPUS_PER_TASK:-0}" -ge 32 ]] || exit 2

root=/public/home/iai806/abacus/agent-runs/20260927-varneb-bto-soft-gridpilot-qx030
code=${root}/code
audit_code=${root}/audit-code-20260927
inputs=${root}/inputs
python=/public/home/iai806/.conda/envs/icu/bin/python
abacus=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus
pseudo_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Pseudopotential
basis_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Orb-DZP-10au
env_script=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh

case "${SLURM_ARRAY_TASK_ID:-}" in
  0) point=q090_q000; preflight_sha=34e2016a4ba0e2ab76a737ecbcd3a0283862247481965382ba2506760e9cebfc ;;
  1) point=q090_q030; preflight_sha=45dea10dc2491d874f24de3d520e7d0de62bc7de62213576bd4410f2f2615aa2 ;;
  *) echo 'unexpected array index' >&2; exit 2 ;;
esac

workdir=${root}/run-${point}
soft_preflight=${root}/preflight-${point}-soft-direction-two-step-20260927.json
comparison=${root}/compare-${point}-curvature-0p05-0p10-27787309.json
runner=${audit_code}/scripts/run_bto_q1q2_soft_probes_abacus.py
output=${root}/soft-probe-${point}-result.json
[[ -s "${soft_preflight}" && -s "${comparison}" && -s "${runner}" \
   && -s "${workdir}/contract.json" && ! -e "${output}" ]] || exit 3
echo "${preflight_sha}  ${soft_preflight}" | sha256sum -c -
echo '4020ab5e4e0d809eb97f4832ba747ff417dff372c465606efc1eb0d512184c68  '"${runner}" | sha256sum -c -

"${python}" - "${comparison}" "${soft_preflight}" <<'PY'
import json, pathlib, sys
comparison, soft = (json.loads(pathlib.Path(p).read_text()) for p in sys.argv[1:])
if not (comparison['status'] == 'two_step_local_curvature_comparison_not_minimum_or_TS_certificate'
        and all(value > 0 for value in comparison['lowest_eigenvalues'])
        and comparison['lowest_direction_comparison']['lowest_mode_absolute_overlap'] > 0.999
        and soft['status'] == 'four_soft_direction_geometries_passed_no_dft'
        and len(soft['probes']) == 4):
    raise SystemExit('audited two-step soft-direction gate failed')
PY

module purge
source "${env_script}"
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS I_MPI_HYDRA_BOOTSTRAP
export PYTHONPATH="${code}" OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
cd "${code}"
"${python}" "${runner}" \
  --preflight "${inputs}/bto_transverse_soft_conditional_${point}_preflight_2026-09-27.json" \
  --soft-preflight "${soft_preflight}" \
  --branch-result "${workdir}/conditional_${point}_result.json" \
  --report "${inputs}/bto_t_to_c_q1q2_strain_gradient_preflight_2026-09-25.json" \
  --reference "${inputs}/cubic_CONTCAR" \
  --force-constants "${inputs}/bto_cubic_gamma_force_constants.npz" \
  --phonopy-eigenpairs "${inputs}/bto_cubic_phonopy_gamma_eigenpairs.npz" \
  --grid-result-manifest "${inputs}/result_manifest_hf_27777454.json" \
  --workdir "${workdir}" --pseudo-dir "${pseudo_dir}" --basis-dir "${basis_dir}" \
  --abacus-bin "${abacus}" --mpi-ranks 32 --output "${output}" --run-dft
