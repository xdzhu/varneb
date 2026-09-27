#!/usr/bin/env bash
# Pre-registered independent BTO (Qz,Qx)=(0.75,0.15) three-start canary.
# Static-only ABACUS; optimization is a separate gate after raw-output audit.
#SBATCH -J varneb-bto-q075-canary
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=00:30:00
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260927-varneb-bto-soft-gridpilot-qx030/slurm-q075-canary-%j.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260927-varneb-bto-soft-gridpilot-qx030/slurm-q075-canary-%j.err

set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 ]] || { echo 'set RUN_DFT=1 explicitly' >&2; exit 2; }
[[ "${SLURM_JOB_PARTITION:-}" == hfacnormal01 && "${SLURM_JOB_NUM_NODES:-0}" == 1 \
   && "${SLURM_NTASKS:-0}" == 1 && "${SLURM_CPUS_PER_TASK:-0}" -ge 32 ]] || exit 2

root=/public/home/iai806/abacus/agent-runs/20260927-varneb-bto-soft-gridpilot-qx030
code=${root}/code
inputs=${root}/inputs
workdir=${root}/run-q075_q015
preflight=${inputs}/bto_transverse_soft_conditional_q075_q015_preflight_2026-09-27.json
python=/public/home/iai806/.conda/envs/icu/bin/python
abacus=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus
pseudo_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Pseudopotential
basis_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Orb-DZP-10au
env_script=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh

[[ -s "${preflight}" && ! -e "${workdir}" ]] || exit 3
echo '1eb33df77355b3c2282333994402110b705318004b5d0305c2fb745967f73ee4  '"${preflight}" | sha256sum -c -
echo 'de0c1c3aeac376b9e4f4ea1f9ff742e85e70f9ba08823824c46c358d547a17bf  '"${code}/examples/run_bto_q1q2_conditional_abacus.py" | sha256sum -c -
echo '91b70d48a488b5a26f1e563bbce1ba0e75581a01d9145158d57108b8dac4bbfc  '"${code}/vcneb/mode_surface.py" | sha256sum -c -

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
  --workdir "${workdir}" --pseudo-dir "${pseudo_dir}" --basis-dir "${basis_dir}" \
  --abacus-bin "${abacus}" --mpi-ranks 32 \
  --gradient-tolerance 0.003 --orthogonal-amplitude-bound 3.0 \
  --max-iterations 2 --max-new-evaluations 3 \
  --result-name canary_starts_result.json --canary-starts-only --run-dft
