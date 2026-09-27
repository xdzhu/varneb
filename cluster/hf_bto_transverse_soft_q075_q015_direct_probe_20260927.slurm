#!/usr/bin/env bash
# Four signed BTO soft-direction statics at fixed Q=(0.75,0.15).
# This cross-check does not certify a conditional PES or a T-to-C barrier.
# ABACUS PBE/100 Ry/Ba-Ti-O 10 au DZP, 4x4x4 electronic k mesh, 32 MPI.
# Gamma force constants originate from the five-atom 1x1x1 cell.
#SBATCH -J varneb-bto-q075-soft4
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=00:45:00
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260927-varneb-bto-soft-gridpilot-qx030/slurm-q075-soft4-%j.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260927-varneb-bto-soft-gridpilot-qx030/slurm-q075-soft4-%j.err

set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 ]] || { echo 'set RUN_DFT=1 explicitly' >&2; exit 2; }
[[ "${SLURM_JOB_PARTITION:-}" == hfacnormal01 && "${SLURM_JOB_NUM_NODES:-0}" == 1 \
   && "${SLURM_NTASKS:-0}" == 1 && "${SLURM_CPUS_PER_TASK:-0}" -ge 32 ]] || exit 2

root=/public/home/iai806/abacus/agent-runs/20260927-varneb-bto-soft-gridpilot-qx030
code=${root}/code
audit_code=${root}/audit-code-20260927
inputs=${root}/inputs
workdir=${root}/run-q075_q015
soft_preflight=${root}/preflight-q075_q015-soft-direction-two-step-20260927.json
runner=${audit_code}/scripts/run_bto_q1q2_soft_probes_abacus.py
output=${root}/soft-probe-q075_q015-result.json
python=/public/home/iai806/.conda/envs/icu/bin/python
abacus=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus
pseudo_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Pseudopotential
basis_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Orb-DZP-10au
env_script=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh

[[ -s "${soft_preflight}" && -s "${runner}" && -s "${workdir}/contract.json" \
   && -s "${workdir}/conditional_q075_q015_result.json" && ! -e "${output}" ]] || exit 3
echo 'da7f47f2cfc47f23cf497257999939d0069dd7a2a30d55ac35bc467ae2b40096  '"${soft_preflight}" | sha256sum -c -
echo '4020ab5e4e0d809eb97f4832ba747ff417dff372c465606efc1eb0d512184c68  '"${runner}" | sha256sum -c -

module purge
source "${env_script}"
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS I_MPI_HYDRA_BOOTSTRAP
export PYTHONPATH="${code}" OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
cd "${code}"
args=(
  --preflight "${inputs}/bto_transverse_soft_conditional_q075_q015_preflight_2026-09-27.json"
  --soft-preflight "${soft_preflight}"
  --branch-result "${workdir}/conditional_q075_q015_result.json"
  --report "${inputs}/bto_t_to_c_q1q2_strain_gradient_preflight_2026-09-25.json"
  --reference "${inputs}/cubic_CONTCAR"
  --force-constants "${inputs}/bto_cubic_gamma_force_constants.npz"
  --phonopy-eigenpairs "${inputs}/bto_cubic_phonopy_gamma_eigenpairs.npz"
  --grid-result-manifest "${inputs}/result_manifest_hf_27777454.json"
  --workdir "${workdir}" --pseudo-dir "${pseudo_dir}" --basis-dir "${basis_dir}"
  --abacus-bin "${abacus}" --mpi-ranks 32 --output "${output}"
)
"${python}" "${runner}" "${args[@]}"
"${python}" "${runner}" "${args[@]}" --run-dft
