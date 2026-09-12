#!/bin/bash
# Single-image ABACUS calculator state-refresh smoke; keep this separate from VCNEB production jobs.
#SBATCH -J abacus-state-refresh
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260912-vcneb-hf/slurm-%x-%j.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260912-vcneb-hf/slurm-%x-%j.err

set -eo pipefail

repo=/public/home/iai806/abacus/agent-runs/20260912-vcneb-hf
env_script=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
python=/public/home/iai806/.conda/envs/icu/bin/python
abacus=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus
pseudo_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Pseudopotential
basis_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/selected_Orbs
workdir=${WORKDIR:-outputs/abacus_state_refresh_smoke}

module purge
source "${env_script}"
export OMP_NUM_THREADS=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi.so
cd "${repo}"

exec "${python}" examples/run_abacus_state_refresh_smoke.py \
  --structure validation/batio3_cubic_to_tetragonal/hf_endpoint_relax_cubic/CONTCAR \
  --workdir "${workdir}" \
  --pseudo-dir "${pseudo_dir}" --basis-dir "${basis_dir}" \
  --command "srun --mpi=pmix_v3 -n ${SLURM_NTASKS} ${abacus}" \
  --ecutwfc "${ECUTWFC:-60}" --scf-thr "${SCF_THR:-1e-6}"
