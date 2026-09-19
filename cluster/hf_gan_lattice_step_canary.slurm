#!/bin/bash
# Full static SCF only: prior image10 reference and independently saved half-step.
#SBATCH -J gan-lattice-step-canary
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
set -eo pipefail
source /public/home/iai806/Software/VASP/env.sh 6.3.2
set -u
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_NTASKS}" == 32 ]] || exit 2
repo=/public/home/iai806/abacus/agent-runs/20260918-varneb-gan-suite
python=/public/home/iai806/.conda/envs/icu/bin/python
root="${repo}/diagnostics/hex_image10_basis/step_canary"
export OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
export VASP_COMMAND="srun --exclusive --exact --nodes=1 --ntasks=32 --ntasks-per-node=32 --cpus-per-task=1 --mpi=pmi2 ${VASP_ROOT}/bin/vasp_std"
cd "${repo}"
"${python}" scripts/probe_vasp_input_contract.py \
  --source-dir cases/b4_b1_hexagonal/input/initial \
  --structures "${root}/inputs/reference.vasp" "${root}/inputs/fractional_step.vasp" \
  --workdir "${root}/full_scf" --symprec 1e-4 --ncores 32 \
  --vasp-bin "${VASP_ROOT}/bin/vasp_std" --timeout 900 --full-scf --run
