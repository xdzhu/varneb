#!/bin/bash
# GaN B4 -> B1 ABACUS/ASE production VC-NEB.
# 29 total images, fixed endpoints, nine concurrent 32-MPI interior workers.
# This is the ABACUS counterpart of the VASP GaN production case.
# Submit only after endpoint identity/static gates have passed.
#SBATCH -J varneb-gan-abacus
#SBATCH -p hfacnormal01
#SBATCH -N 3
#SBATCH --ntasks=291
#SBATCH --ntasks-per-node=97
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260921-varneb-material/slurm-%x-%j.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260921-varneb-material/slurm-%x-%j.err

set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 ]] || { echo "set RUN_DFT=1" >&2; exit 2; }

repo=${REPO:-/public/home/iai806/abacus/agent-runs/20260921-varneb-material}
python=${PYTHON:-/public/home/iai806/.conda/envs/icu/bin/python}
abacus=${ABACUS_BIN:-/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus}
env_script=${ABACUS_ENV:-/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh}
pseudo_dir=${PSEUDO_DIR:-/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Pseudopotential}
basis_dir=${BASIS_DIR:-/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Orb-DZP-10au}
initial=${INITIAL:-cases/gan/endpoints/B4_POSCAR}
final=${FINAL:-cases/gan/endpoints/B1_POSCAR}
workdir=${WORKDIR:-cases/gan/abacus_vcneb_production}
n_images=${N_IMAGES:-29}
image_workers=${IMAGE_WORKERS:-9}
image_mpi=${IMAGE_MPI:-32}
steps=${STEPS:-300}
fmax=${FMAX:-0.10}

(( n_images == 29 && image_workers == 9 && image_mpi == 32 )) || {
  echo "GaN production requires 29 images, 9 concurrent workers, and 32 MPI ranks per worker" >&2
  exit 2
}
(( SLURM_NTASKS >= image_workers * image_mpi )) || exit 2

module purge
source "${env_script}"
export OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
cd "${repo}"

exec "${python}" examples/run_vcneb_abacus.py \
  --initial "${initial}" --final "${final}" --workdir "${workdir}" \
  --n-images "${n_images}" --fmax "${fmax}" --steps "${steps}" --k 0.20 \
  --optimizer FIRE --image-workers "${image_workers}" --no-climb \
  --command "srun --exclusive --mpi=pmix_v3 -n ${image_mpi} ${abacus}" \
  --pseudo-dir "${pseudo_dir}" --basis-dir "${basis_dir}" \
  --pp Ga=Ga.upf --pp N=N.upf \
  --basis Ga=Ga_gga_10au_100Ry_2s2p2d1f.orb \
  --basis N=N_gga_10au_100Ry_2s2p1d.orb \
  --ecutwfc 100 --scf-thr 1e-8 --scf-nmax 150 --mixing-beta 0.3 \
  --kpts 4 4 3 --mic --cell-interpolation log_strain --mapping auto \
  --align-translation --minimum-distance 1.4 --maximum-deformation 0.50
