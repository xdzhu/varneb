#!/bin/bash
# BTO T->C VASP VCNEB: one manager plus five independent 32-MPI interiors.
# This template never guesses a licensed VASP installation or PAW dataset.
# Set RUN_DFT=1, VASP_BIN, VASP_INITIAL_DIR and VASP_FINAL_DIR explicitly.
# N_IMAGES includes endpoints, so seven total images means five workers.
#SBATCH -J varneb-bto-vasp
#SBATCH -p hfacnormal01
#SBATCH -N 2
#SBATCH --ntasks=160
#SBATCH --cpus-per-task=1
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260912-vcneb-hf/slurm-%x-%j.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260912-vcneb-hf/slurm-%x-%j.err

set -euo pipefail
source /public/home/iai806/Software/VASP/env.sh 6.3.2

repo=${REPO:-/public/home/iai806/abacus/agent-runs/20260912-vcneb-hf}
python=${PYTHON:-/public/home/iai806/.conda/envs/icu/bin/python}
n_images=${N_IMAGES:-7}
n_interior=$((n_images - 2))
image_workers=${IMAGE_WORKERS:-${n_interior}}
image_mpi=${IMAGE_MPI:-32}
steps=${STEPS:-300}
fmax=${FMAX:-0.10}
spring=${SPRING:-0.2}
optimizer=${OPTIMIZER:-FIRE}
initial_dir=${VASP_INITIAL_DIR:-}
final_dir=${VASP_FINAL_DIR:-}
endpoint_identity_gate=${ENDPOINT_IDENTITY_GATE:-}
static_baseline=${VASP_STATIC_BASELINE:-}
workdir=${WORKDIR:-validation/batio3_cubic_to_tetragonal/vcneb_tetragonal_to_cubic_n${n_images}_${optimizer,,}_vasp}

if [[ "${RUN_DFT:-0}" != 1 ]]; then
  echo "Refusing a 160-rank allocation without RUN_DFT=1; use the local --validate-only driver first" >&2
  exit 2
fi
if (( n_images != 7 || image_workers != n_interior || image_workers * image_mpi != 160 )); then
  echo "BTO VASP validation is fixed at 7 total images and five 32-MPI interior workers" >&2
  exit 2
fi
if (( SLURM_NTASKS < image_workers * image_mpi )); then
  echo "allocation is insufficient for five 32-MPI VASP workers" >&2
  exit 2
fi
: "${VASP_BIN:?Set VASP_BIN to the licensed VASP executable visible on hfacnormal01}"
: "${VASP_INITIAL_DIR:?Set VASP_INITIAL_DIR (CONTCAR, INCAR, KPOINTS, POTCAR)}"
: "${VASP_FINAL_DIR:?Set VASP_FINAL_DIR (CONTCAR)}"
: "${ENDPOINT_IDENTITY_GATE:?Set ENDPOINT_IDENTITY_GATE to a passing endpoint_identity_gate.json}"
: "${VASP_STATIC_BASELINE:?Set VASP_STATIC_BASELINE to a completed vasp_static_summary.json}"
for path in "${initial_dir}/CONTCAR" "${initial_dir}/INCAR" "${initial_dir}/KPOINTS" "${initial_dir}/POTCAR" "${final_dir}/CONTCAR"; do
  [[ -f "${path}" ]] || { echo "missing required VASP input: ${path}" >&2; exit 2; }
done

export OMP_NUM_THREADS=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
export VCNEB_GIT_REVISION="${VCNEB_GIT_REVISION:-remote-sync-unknown}"
export VASP_COMMAND="srun --exclusive --exact --nodes=1 --ntasks=${image_mpi} --ntasks-per-node=${image_mpi} --cpus-per-task=1 --distribution=block:block --mpi=pmi2 ${VASP_BIN}"
cd "${repo}"
"${python}" scripts/validate_endpoint_identity_gate.py --gate "${endpoint_identity_gate}"
"${python}" scripts/validate_vasp_static_baseline_gate.py \
  --summary "${static_baseline}" --endpoint-gate "${endpoint_identity_gate}" \
  --initial-dir "${initial_dir}" --git-revision "${VCNEB_GIT_REVISION}"

exec "${python}" examples/run_vcneb_vasp.py \
  --initial "${initial_dir}" --final "${final_dir}" --workdir "${workdir}" \
  --n-images "${n_images}" --fmax "${fmax}" --steps "${steps}" --k "${spring}" \
  --optimizer "${optimizer}" --image-workers "${image_workers}" \
  --no-climb --mic --cell-interpolation log_strain --mapping auto --align-translation \
  --minimum-distance 1.6 --maximum-deformation 0.10
