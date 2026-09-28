#!/bin/bash
# Resume the 22-point, independently audited Qy=0 cache after its optimizer
# coordinate amplitude (not a DFT parameter) hit the 3.0 upper bound.
#SBATCH -J varneb-bto-qy0-Tedge
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=05:00:00
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260928-varneb-bto-qy0-pilot/slurm-q120-bound-%j.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260928-varneb-bto-qy0-pilot/slurm-q120-bound-%j.err

set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 ]] || { echo 'set RUN_DFT=1 explicitly' >&2; exit 2; }
[[ "${SLURM_JOB_PARTITION:-}" == hfacnormal01 && "${SLURM_JOB_NUM_NODES:-0}" == 1 \
   && "${SLURM_NTASKS:-0}" == 1 && "${SLURM_CPUS_PER_TASK:-0}" -ge 32 ]] || {
  echo 'requires one 32-CPU Slurm task on hfacnormal01' >&2
  exit 2
}

root=/public/home/iai806/abacus/agent-runs/20260928-varneb-bto-qy0-pilot
code=${root}/code-bound-continuation
inputs=${root}/inputs
workdir=${root}/run-q120_q030
audit=${root}/audit-q120_q030-27794649_3.json
preflight=${inputs}/bto_qy0_q120_q030_preflight_2026-09-28.json
result=q120_q030_restricted_bound4_result.json
python=/public/home/iai806/.conda/envs/icu/bin/python
abacus=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus
pseudo_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Pseudopotential
basis_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Orb-DZP-10au
env_script=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh

[[ -s "${workdir}/contract.json" && -s "${audit}" && ! -e "${workdir}/${result}" ]] || {
  echo 'audited source cache required; continuation result must be unused' >&2
  exit 2
}
echo '29d68588c80959106557d7672ee8008519541891541a242a5b06a63b14c1b3e6  '${preflight} | sha256sum -c -
echo 'd92f7c9de77f9c78f691092b40f94c66d6d3b9a6ccddcfdd718d8d644b76e722  '${audit} | sha256sum -c -
echo '42b9e1d5ebecf096261a85c79606ce62843e40a92e51e286821f74f44a148258  '${code}/examples/run_bto_q1q2_conditional_abacus.py | sha256sum -c -
echo 'bffce25d9054494cc413e8e530d7e47248a847aa56e05da91ce3f11d95d59107  '${code}/scripts/audit_bto_qy0_failed_cache.py | sha256sum -c -
echo '55bcf70c3af12b4a5d4eb445ae30f8461c07f18890afc34c05ee3b084b64164a  '${code}/vcneb/mode_surface.py | sha256sum -c -

module purge
source "${env_script}"
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS I_MPI_HYDRA_BOOTSTRAP
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
cd "${code}"
"${python}" examples/run_bto_q1q2_conditional_abacus.py \
  --preflight "${preflight}" \
  --report "${inputs}/bto_t_to_c_q1q2_strain_gradient_preflight_2026-09-25.json" \
  --reference "${inputs}/cubic_CONTCAR" \
  --force-constants "${inputs}/bto_cubic_gamma_force_constants.npz" \
  --phonopy-eigenpairs "${inputs}/bto_cubic_phonopy_gamma_eigenpairs.npz" \
  --gamma-provenance "${inputs}/bto_cubic_gamma_phonon_provenance.json" \
  --force-sets "${inputs}/bto_cubic_gamma_FORCE_SETS" \
  --eigenpairs-provenance "${inputs}/bto_cubic_phonopy_gamma_eigenpairs_provenance.json" \
  --grid-result-manifest "${inputs}/result_manifest_hf_27777454.json" \
  --workdir "${workdir}" \
  --pseudo-dir "${pseudo_dir}" --basis-dir "${basis_dir}" \
  --abacus-bin "${abacus}" --mpi-ranks 32 \
  --gradient-tolerance 0.003 --orthogonal-amplitude-bound 4.0 \
  --max-iterations 64 --max-new-evaluations 120 \
  --max-absolute-stress-kbar 2.0 \
  --start-from-failed-cache-audit "${audit}" \
  --result-name "${result}" --resume --run-dft
