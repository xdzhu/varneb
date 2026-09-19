#!/bin/bash
#SBATCH -J varneb-gan-endpoint
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
set -eo pipefail
source /public/home/iai806/Software/VASP/env.sh 6.3.2
set -u
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_NTASKS}" == 32 ]] || exit 2
[[ "${ENDPOINT}" == b3 || "${ENDPOINT}" == b1 ]] || exit 2
repo=/public/home/iai806/abacus/agent-runs/20260918-varneb-gan-suite
python=/public/home/iai806/.conda/envs/icu/bin/python
export OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
export VASP_COMMAND="srun --exclusive --exact --nodes=1 --ntasks=32 --ntasks-per-node=32 --cpus-per-task=1 --mpi=pmi2 ${VASP_ROOT}/bin/vasp_std"
cd "${repo}"
"${python}" scripts/relax_vasp_endpoint.py \
  --structure "cases/endpoint_seeds/${ENDPOINT}/POSCAR" \
  --source-dir "cases/endpoint_seeds/${ENDPOINT}" \
  --workdir "cases/endpoint_relax/${ENDPOINT}" \
  --pressure-gpa 45.0 --fmax 0.02 --steps 120 --maxstep 0.05 \
  --symprec 1e-4 --ncores 32 --vasp-bin "${VASP_ROOT}/bin/vasp_std"
