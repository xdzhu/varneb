#!/bin/bash
# Distributed VCNEB manager plus one concurrent ABACUS worker per interior image.
#
# Seven total HfO2 images have five interior images, so they use five 32-MPI
# workers (160 ranks) spread across two 128-CPU hfacnormal01 nodes.  Endpoints
# are evaluated once by the Python manager and cached; only interiors launch
# isolated Slurm job steps on each VC-NEB iteration.
#SBATCH -J vcneb-hfo2-dist
#SBATCH -p hfacnormal01
#SBATCH -N 2
#SBATCH --ntasks=160
# Five 32-rank workers need at most three workers on one node and two on the
# other; reserve 96 task slots per node instead of the full 128-core nodes.
#SBATCH --ntasks-per-node=96
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
n_interior=$((n_images - 2))
image_workers=${IMAGE_WORKERS:-${n_interior}}
image_mpi=${IMAGE_MPI:-32}
node_capacity=${TASKS_PER_NODE:-128}
image_retries=${IMAGE_RETRIES:-0}
image_cache_dir=${IMAGE_CACHE_DIR:-}
image_cache_namespace=${IMAGE_CACHE_NAMESPACE:-}
steps=${STEPS:-300}
fmax=${FMAX:-0.05}
spring=${SPRING:-0.2}
optimizer=${OPTIMIZER:-FIRE}
maxstep=${MAXSTEP:-0.02}
climb_after=${CLIMB_AFTER:-100}
cell_interpolation=${CELL_INTERPOLATION:-log_strain}
mapping=${MAPPING:-auto}
align_translation=${ALIGN_TRANSLATION:-1}
ecutwfc=${ECUTWFC:-100}
scf_thr=${SCF_THR:-1e-8}
scf_nmax=${SCF_NMAX:-150}
mixing_beta=${MIXING_BETA:-0.3}
kpt_x=${KPT_X:-2}
kpt_y=${KPT_Y:-2}
kpt_z=${KPT_Z:-2}
workdir=${WORKDIR:-validation/hfo2_t_to_po/vcneb_n${n_images}_${optimizer,,}_distributed}
initial_endpoint=${INITIAL_ENDPOINT:-validation/hfo2_t_to_po/relaxed_T/CONTCAR}
final_endpoint=${FINAL_ENDPOINT:-validation/hfo2_t_to_po/relaxed_PO/CONTCAR}

if (( n_images < 3 )); then
  echo "N_IMAGES must be at least 3 (two fixed endpoints plus an interior)" >&2
  exit 2
fi
if (( image_workers < 1 || image_mpi < 1 || node_capacity < image_mpi )); then
  echo "IMAGE_WORKERS, IMAGE_MPI and TASKS_PER_NODE are invalid" >&2
  exit 2
fi
if (( image_workers > n_interior )); then
  echo "IMAGE_WORKERS=${image_workers} cannot exceed interior images=${n_interior}" >&2
  exit 2
fi
required_tasks=$((image_workers * image_mpi))
required_nodes=$(( (required_tasks + node_capacity - 1) / node_capacity ))
if (( required_tasks > SLURM_NTASKS )); then
  echo "worker ranks exceed allocation: ${required_tasks} > ${SLURM_NTASKS}" >&2
  exit 2
fi
if (( required_nodes > SLURM_JOB_NUM_NODES )); then
  echo "worker nodes exceed allocation: ${required_nodes} > ${SLURM_JOB_NUM_NODES}" >&2
  exit 2
fi
case "${cell_interpolation}" in
  linear|log_strain) ;;
  *) echo "CELL_INTERPOLATION must be linear or log_strain" >&2; exit 2 ;;
esac
case "${mapping}" in
  identity|auto) ;;
  *) echo "MAPPING must be identity or auto" >&2; exit 2 ;;
esac

module purge
source "${env_script}"
export OMP_NUM_THREADS=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi.so
cd "${repo}"
export VCNEB_GIT_REVISION="${VCNEB_GIT_REVISION:-remote-sync-unknown}"

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

cache_args=()
if [[ -n "${image_cache_dir}" ]]; then
  cache_args+=(--image-cache-dir "${image_cache_dir}")
  if [[ -n "${image_cache_namespace}" ]]; then
    cache_args+=(--image-cache-namespace "${image_cache_namespace}")
  fi
fi

path_args=(--mic --cell-interpolation "${cell_interpolation}" --mapping "${mapping}")
if [[ "${align_translation}" != 0 ]]; then
  path_args+=(--align-translation)
fi

# --nodes=1 is essential: without it Slurm may spread a 32-rank image worker
# over both nodes, destroying the intended interior-image concurrency.
worker_command="srun --exclusive --nodes=1 --ntasks=${image_mpi} --ntasks-per-node=${image_mpi} --cpus-per-task=1 --distribution=block:block --mpi=pmix_v3 ${abacus}"

exec "${python}" examples/run_vcneb_abacus.py \
  --initial "${initial_endpoint}" \
  --final "${final_endpoint}" \
  --workdir "${workdir}" --n-images "${n_images}" --fmax "${fmax}" \
  --steps "${steps}" --k "${spring}" --optimizer "${optimizer}" \
  --maxstep "${maxstep}" --image-workers "${image_workers}" \
  --image-retries "${image_retries}" \
  --image-manifest "${workdir}/image_worker_manifest.jsonl" \
  --command "${worker_command}" \
  --pseudo-dir "${pseudo_dir}" --basis-dir "${basis_dir}" \
  --pp Hf=Hf.upf --pp O=O.upf \
  --basis Hf=Hf_gga_10au_100Ry_4s2p2d1f.orb \
  --basis O=O_gga_10au_100Ry_2s2p1d.orb \
  --ecutwfc "${ecutwfc}" --scf-thr "${scf_thr}" --scf-nmax "${scf_nmax}" \
  --mixing-beta "${mixing_beta}" --kpts "${kpt_x}" "${kpt_y}" "${kpt_z}" \
  "${path_args[@]}" \
  --minimum-distance 1.6 --maximum-deformation 0.25 \
  --climb-after "${climb_after}" \
  "${cache_args[@]}" \
  "${climb_args[@]}" "${resume_args[@]}"
