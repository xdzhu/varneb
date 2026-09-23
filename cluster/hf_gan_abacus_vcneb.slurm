#!/bin/bash
# GaN B4 -> B1 ABACUS/ASE production VC-NEB.
# 29 total images (27 interior), fixed endpoints, nine concurrent 4-MPI workers.
# Intel MPI's mpirun starts one Hydra/Slurm proxy per worker.  Give each proxy
# four CPUs; a one-CPU Slurm task pins all four ABACUS ranks to one core.
# This is the ABACUS counterpart of the VASP GaN production case.
# Submit only after endpoint identity/static gates have passed.
#SBATCH -J varneb-gan-abacus
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=10
#SBATCH --ntasks-per-node=10
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260921-varneb-material/slurm-%x-%j.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260921-varneb-material/slurm-%x-%j.err

set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 ]] || { echo "set RUN_DFT=1" >&2; exit 2; }

repo=${REPO:-/public/home/iai806/abacus/agent-runs/20260921-varneb-material/varneb_abacus_45p7_source}
python=${PYTHON:-/public/home/iai806/.conda/envs/icu/bin/python}
abacus=${ABACUS_BIN:-/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus}
env_script=${ABACUS_ENV:-/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh}
pseudo_dir=${PSEUDO_DIR:-/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Pseudopotential}
basis_dir=${BASIS_DIR:-/public/home/iai806/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Orb-DZP-10au}
initial=${INITIAL:?set INITIAL to the pressure-relaxed B4 endpoint}
final=${FINAL:?set FINAL to the pressure-relaxed B1 endpoint}
workdir=${WORKDIR:?set WORKDIR to a new production directory}
pressure_gpa=${PRESSURE_GPA:-45.7}
n_images=${N_IMAGES:-29}
image_workers=${IMAGE_WORKERS:-9}
image_mpi=${IMAGE_MPI:-4}
steps=${STEPS:-300}
fmax=${FMAX:-0.10}
maximum_cell_step=${MAXIMUM_CELL_STEP:-0.05}
candidate_step_retries=${CANDIDATE_STEP_RETRIES:-8}
resume=${RESUME:-0}
resume_trajectory=${RESUME_TRAJECTORY:-}
resume_step=${RESUME_STEP:-0}
endpoint_summary=${ENDPOINT_STATIC_SUMMARY:-}

(( n_images == 29 && image_workers == 9 && image_mpi == 4 )) || {
  echo "GaN production requires 29 images, 9 concurrent workers, and 4 MPI ranks per worker" >&2
  exit 2
}
[[ "${pressure_gpa}" == 45.7 ]] || {
  echo "this GaN B4->B1 case is defined at 45.7 GPa" >&2
  exit 2
}
(( SLURM_JOB_NUM_NODES == 1 && SLURM_NTASKS >= image_workers + 1 \
   && SLURM_CPUS_PER_TASK >= image_mpi )) || {
  echo "need one node with one manager task plus nine 4-CPU MPI bootstrap tasks" >&2
  exit 2
}

module purge
source "${env_script}"
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS I_MPI_HYDRA_BOOTSTRAP
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONUNBUFFERED=1
cd "${repo}"
[[ -s "${initial}" && -s "${final}" ]] || {
  echo "pressure-relaxed endpoint files are missing" >&2; exit 2;
}
if [[ "${resume}" != 1 ]]; then
  [[ ! -e "${workdir}" ]] || { echo "refusing existing WORKDIR=${workdir}" >&2; exit 2; }
fi
mkdir -p "${workdir}"

if [[ "${REQUIRE_ENDPOINT_STATIC_GATE:-1}" == 1 ]]; then
  [[ -n "${endpoint_summary}" ]] || {
    echo "ENDPOINT_STATIC_SUMMARY is required before production ABACUS VC-NEB" >&2
    exit 2
  }
  if [[ "${endpoint_summary}" != /* ]]; then
    endpoint_summary="${repo}/${endpoint_summary}"
  fi
  "${python}" scripts/validate_ase_static_gate.py "${endpoint_summary}" \
    --fmax "${ENDPOINT_FMAX:-0.10}" \
    --stress-kbar "${ENDPOINT_STRESS_KBAR:-2.0}" \
    --pressure-gpa "${pressure_gpa}" \
    --output "${workdir}/endpoint_static_gate.json"
fi

command=("${python}" examples/run_vcneb_abacus.py \
  --initial "${initial}" --final "${final}" --workdir "${workdir}" \
  --n-images "${n_images}" --fmax "${fmax}" --steps "${steps}" --k 0.20 \
  --pressure-gpa "${pressure_gpa}" \
  --optimizer FIRE --image-workers "${image_workers}" --no-climb \
  --candidate-step-retries "${candidate_step_retries}" --maximum-cell-step "${maximum_cell_step}" \
  --command "mpirun -np ${image_mpi} ${abacus}" \
  --pseudo-dir "${pseudo_dir}" --basis-dir "${basis_dir}" \
  --pp Ga=Ga.upf --pp N=N.upf \
  --basis Ga=Ga_gga_10au_100Ry_2s2p2d1f.orb \
  --basis N=N_gga_10au_100Ry_2s2p1d.orb \
  --ecutwfc 100 --scf-thr 1e-8 --scf-nmax 150 --mixing-beta 0.3 \
  --kpts 4 4 3 --mic --cell-interpolation log_strain --mapping auto \
  --align-translation --minimum-distance 1.4 --maximum-deformation 0.50)
if [[ -n "${endpoint_summary}" ]]; then
  command+=(--endpoint-static-summary "${endpoint_summary}")
fi
if [[ "${resume}" == 1 ]]; then
  command+=(--resume --resume-step "${resume_step}")
  if [[ -n "${resume_trajectory}" ]]; then
    command+=(--resume-trajectory "${resume_trajectory}")
  fi
fi
exec "${command[@]}"
