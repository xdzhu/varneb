#!/bin/bash
# One VCNEB controller with concurrent image-level Slurm job steps.
# Four 32-MPI image workers fill one 128-CPU hfacnormal01 node.
#SBATCH -J vcneb-bto-par
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260912-vcneb-hf/slurm-%x-%j.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260912-vcneb-hf/slurm-%x-%j.err

set -eo pipefail

repo=/public/home/iai806/abacus/agent-runs/20260912-vcneb-hf
env_script=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
python=/public/home/iai806/.conda/envs/icu/bin/python
abacus=/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus
pseudo_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Pseudopotential
basis_dir=/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Orb-DZP-10au

n_images=${N_IMAGES:-7}
image_workers=${IMAGE_WORKERS:-4}
image_mpi=${IMAGE_MPI:-32}
image_retries=${IMAGE_RETRIES:-0}
steps=${STEPS:-300}
fmax=${FMAX:-0.03}
spring=${SPRING:-0.2}
optimizer=${OPTIMIZER:-FIRE}
maxstep=${MAXSTEP:-0.02}
climb_after=${CLIMB_AFTER:-100}
ecutwfc=${ECUTWFC:-100}
scf_thr=${SCF_THR:-1e-8}
scf_nmax=${SCF_NMAX:-150}
mixing_beta=${MIXING_BETA:-0.3}
kpt_x=${KPT_X:-4}
kpt_y=${KPT_Y:-4}
kpt_z=${KPT_Z:-4}
direction=${DIRECTION:-tetragonal_to_cubic}
endpoint_tag=${ENDPOINT_TAG:-pbe100_dzp10au}
workdir=${WORKDIR:-validation/batio3_cubic_to_tetragonal/vcneb_${direction}_n${n_images}_${optimizer,,}_${endpoint_tag}_parallel}

if (( image_workers < 1 || image_mpi < 1 )); then
  echo "IMAGE_WORKERS and IMAGE_MPI must be positive" >&2
  exit 2
fi
if (( image_workers * image_mpi > SLURM_NTASKS )); then
  echo "IMAGE_WORKERS*IMAGE_MPI exceeds allocated tasks: ${image_workers}*${image_mpi} > ${SLURM_NTASKS}" >&2
  exit 2
fi

case "${direction}" in
  tetragonal_to_cubic)
    initial_endpoint="tetragonal"
    final_endpoint="cubic"
    ;;
  cubic_to_tetragonal)
    initial_endpoint="cubic"
    final_endpoint="tetragonal"
    ;;
  *)
    echo "unsupported DIRECTION=${direction}" >&2
    exit 2
    ;;
esac

module purge
source "${env_script}"
export OMP_NUM_THREADS=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi.so
cd "${repo}"

climb_args=()
if [[ "${NO_CLIMB:-0}" != 0 ]]; then
  climb_args+=(--no-climb)
fi
resume_args=()
if [[ "${RESUME:-0}" != 0 ]]; then
  resume_args+=(--resume)
  if [[ -n "${RESUME_TRAJECTORY:-}" ]]; then
    resume_args+=(--resume-trajectory "${RESUME_TRAJECTORY}")
  fi
  if [[ -n "${RESUME_STEP:-}" ]]; then
    resume_args+=(--resume-step "${RESUME_STEP}")
  fi
fi

exec "${python}" examples/run_vcneb_abacus.py \
  --initial "validation/batio3_cubic_to_tetragonal/endpoint_relax_${initial_endpoint}_${endpoint_tag}/CONTCAR" \
  --final "validation/batio3_cubic_to_tetragonal/endpoint_relax_${final_endpoint}_${endpoint_tag}/CONTCAR" \
  --workdir "${workdir}" --n-images "${n_images}" --fmax "${fmax}" \
  --steps "${steps}" --k "${spring}" --optimizer "${optimizer}" \
  --maxstep "${maxstep}" --image-workers "${image_workers}" --image-retries "${image_retries}" \
  --command "srun --exclusive --mpi=pmix_v3 -n ${image_mpi} ${abacus}" \
  --pseudo-dir "${pseudo_dir}" --basis-dir "${basis_dir}" \
  --pp Ba=Ba.upf --pp Ti=Ti.upf --pp O=O.upf \
  --basis Ba=Ba_gga_10au_100Ry_4s2p1d.orb \
  --basis Ti=Ti_gga_10au_100Ry_4s2p2d1f.orb \
  --basis O=O_gga_10au_100Ry_2s2p1d.orb \
  --ecutwfc "${ecutwfc}" --scf-thr "${scf_thr}" --scf-nmax "${scf_nmax}" \
  --mixing-beta "${mixing_beta}" --kpts "${kpt_x}" "${kpt_y}" "${kpt_z}" \
  --mic --cell-interpolation log_strain --mapping auto --align-translation \
  --minimum-distance 1.6 --maximum-deformation 0.10 \
  --climb-after "${climb_after}" \
  "${climb_args[@]}" "${resume_args[@]}"
