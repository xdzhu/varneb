#!/bin/bash
# Real fixed-endpoint VASP SCF baseline for the BTO cross-calculator case.
# It consumes one 32-MPI allocation and evaluates only endpoint 00. It is not
# an NEB step and never lets VASP optimize ions or the cell.
#SBATCH -J varneb-bto-vasp-static
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260912-vcneb-hf/slurm-%x-%j.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260912-vcneb-hf/slurm-%x-%j.err

set -euo pipefail

repo=${REPO:-/public/home/iai806/abacus/agent-runs/20260912-vcneb-hf}
python=${PYTHON:-/public/home/iai806/.conda/envs/icu/bin/python}
initial_dir=${VASP_INITIAL_DIR:-}
final_dir=${VASP_FINAL_DIR:-}
endpoint_identity_gate=${ENDPOINT_IDENTITY_GATE:-}
workdir=${WORKDIR:-validation/batio3_cubic_to_tetragonal/vasp_static_baseline}

if [[ "${RUN_DFT:-0}" != 1 ]]; then
  echo "Refusing a 32-rank allocation without RUN_DFT=1" >&2
  exit 2
fi
: "${VASP_BIN:?Set VASP_BIN to the licensed VASP executable visible on hfacnormal01}"
: "${VASP_INITIAL_DIR:?Set VASP_INITIAL_DIR (CONTCAR, INCAR, KPOINTS, POTCAR)}"
: "${VASP_FINAL_DIR:?Set VASP_FINAL_DIR (CONTCAR)}"
: "${ENDPOINT_IDENTITY_GATE:?Set ENDPOINT_IDENTITY_GATE to a passing endpoint_identity_gate.json}"
for path in "${initial_dir}/CONTCAR" "${initial_dir}/INCAR" "${initial_dir}/KPOINTS" "${initial_dir}/POTCAR" "${final_dir}/CONTCAR"; do
  [[ -f "${path}" ]] || { echo "missing required VASP input: ${path}" >&2; exit 2; }
done

export OMP_NUM_THREADS=1
export VCNEB_GIT_REVISION="${VCNEB_GIT_REVISION:-remote-sync-unknown}"
export VASP_COMMAND="srun --exclusive --nodes=1 --ntasks=32 --ntasks-per-node=32 ${VASP_BIN}"
cd "${repo}"
"${python}" scripts/validate_endpoint_identity_gate.py --gate "${endpoint_identity_gate}"

exec "${python}" examples/run_vcneb_vasp.py \
  --initial "${initial_dir}" --final "${final_dir}" --workdir "${workdir}" \
  --n-images 7 --fmax 0.10 --steps 0 --k 0.2 --image-workers 0 \
  --no-climb --mic --cell-interpolation log_strain --mapping auto --align-translation \
  --minimum-distance 1.6 --maximum-deformation 0.10 --static-only
