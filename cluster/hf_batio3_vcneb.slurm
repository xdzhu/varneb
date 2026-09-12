#!/bin/bash
#SBATCH -J vcneb-bto-neb
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=16
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

n_images=${N_IMAGES:-7}
steps=${STEPS:-300}
fmax=${FMAX:-0.03}
spring=${SPRING:-0.2}
optimizer=${OPTIMIZER:-FIRE}
maxstep=${MAXSTEP:-0.02}
workdir=${WORKDIR:-validation/batio3_cubic_to_tetragonal/vcneb_n${n_images}_${optimizer,,}_pbe100}

module purge
source "${env_script}"
export OMP_NUM_THREADS=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi.so
cd "${repo}"

climb_args=()
if [[ "${NO_CLIMB:-1}" != 0 ]]; then
  climb_args+=(--no-climb)
fi
resume_args=()
if [[ "${RESUME:-0}" != 0 ]]; then
  resume_args+=(--resume)
fi

exec "${python}" examples/run_vcneb_abacus.py \
  --initial validation/batio3_cubic_to_tetragonal/endpoint_relax_cubic/CONTCAR \
  --final validation/batio3_cubic_to_tetragonal/endpoint_relax_tetragonal/CONTCAR \
  --workdir "${workdir}" --n-images "${n_images}" --fmax "${fmax}" \
  --steps "${steps}" --k "${spring}" --optimizer "${optimizer}" \
  --maxstep "${maxstep}" --command "srun --mpi=pmix_v3 -n ${SLURM_NTASKS} ${abacus}" \
  --pseudo-dir "${pseudo_dir}" --basis-dir "${basis_dir}" \
  --pp Ba=Ba.upf --pp Ti=Ti.upf --pp O=O.upf \
  --basis Ba=Ba_gga_10au_100Ry_4s2p1d.orb \
  --basis Ti=Ti_gga_8au_100Ry_4s2p2d1f.orb \
  --basis O=O_gga_7au_100Ry_2s2p1d.orb \
  --ecutwfc 100 --scf-thr 1e-8 --scf-nmax 150 --mixing-beta 0.3 \
  --kpts 4 4 4 --mic --cell-interpolation log_strain --mapping auto --align-translation \
  "${climb_args[@]}" "${resume_args[@]}"
