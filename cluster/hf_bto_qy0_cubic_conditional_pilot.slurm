#!/bin/bash
# One fresh BTO cubic-point canary on the Q_y=0 symmetry-restricted sheet.
# This is a local branch/endpoint check, not a complete conditional PES.
#SBATCH -J varneb-bto-qy0-c
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=02:00:00
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260928-varneb-bto-qy0-pilot/slurm-%j.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260928-varneb-bto-qy0-pilot/slurm-%j.err

set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 ]] || { echo 'set RUN_DFT=1 explicitly' >&2; exit 2; }
[[ "${SLURM_JOB_PARTITION:-}" == hfacnormal01 && "${SLURM_JOB_NUM_NODES:-0}" == 1 \
   && "${SLURM_NTASKS:-0}" == 1 && "${SLURM_CPUS_PER_TASK:-0}" -ge 32 ]] || {
  echo 'requires one 32-CPU Slurm task on hfacnormal01' >&2
  exit 2
}

root=/public/home/iai806/abacus/agent-runs/20260928-varneb-bto-qy0-pilot
code=${root}/code
inputs=${root}/inputs
preflight=${inputs}/bto_qy0_q000_q000_preflight_2026-09-28.json
python=/public/home/iai806/.conda/envs/icu/bin/python
abacus=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus
pseudo_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Pseudopotential
basis_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Orb-DZP-10au
env_script=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh

[[ -d "${code}" && -d "${inputs}" && ! -e "${root}/run-q000" ]] || {
  echo 'staged inputs/code are required and the run directory must be unused' >&2
  exit 2
}
echo "bca070dea043423966cae79bbecbfe97aeb2474fec1a6d6b45f6320f41be5a73  ${preflight}" | sha256sum -c -
echo "4db5823766b88c4e9e103827bee198e51a9a44954bfb532ade6ff5cdcf6c3f1a  ${code}/examples/run_bto_q1q2_conditional_abacus.py" | sha256sum -c -
echo "55bcf70c3af12b4a5d4eb445ae30f8461c07f18890afc34c05ee3b084b64164a  ${code}/vcneb/mode_surface.py" | sha256sum -c -

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
  --workdir "${root}/run-q000" \
  --pseudo-dir "${pseudo_dir}" --basis-dir "${basis_dir}" \
  --abacus-bin "${abacus}" --mpi-ranks 32 \
  --gradient-tolerance 0.003 --orthogonal-amplitude-bound 3.0 \
  --max-iterations 30 --max-new-evaluations 60 \
  --max-absolute-stress-kbar 2.0 \
  --result-name q000_restricted_result.json --run-dft
