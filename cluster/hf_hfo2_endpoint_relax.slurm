#!/bin/bash
#SBATCH -J hfo2-endpoint
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
structure=${STRUCTURE:?STRUCTURE must be set}
workdir=${WORKDIR:?WORKDIR must be set}

module purge
source "${env_script}"
export OMP_NUM_THREADS=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi.so
cd "${repo}"

exec "${python}" examples/relax_abacus_endpoint.py \
  --structure "${structure}" --workdir "${workdir}" \
  --pseudo-dir "${pseudo_dir}" --basis-dir "${basis_dir}" \
  --command "srun --mpi=pmix_v3 -n ${SLURM_NTASKS} ${abacus}" \
  --pp Hf=Hf.upf --pp O=O.upf \
  --basis Hf=Hf_gga_7au_100Ry_4s2p2d1f.orb \
  --basis O=O_gga_7au_100Ry_2s2p1d.orb \
  --ecutwfc "${ECUTWFC:-60}" --scf-thr "${SCF_THR:-1e-6}" \
  --scf-nmax "${SCF_NMAX:-100}" --mixing-beta "${MIXING_BETA:-0.7}" \
  --kpts "${KPT_X:-1}" "${KPT_Y:-1}" "${KPT_Z:-1}" \
  --fmax "${FMAX:-0.05}" --steps "${STEPS:-100}" \
  --optimizer "${OPTIMIZER:-FIRE}" --maxstep "${MAXSTEP:-0.02}"
