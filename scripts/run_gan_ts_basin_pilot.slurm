#!/bin/bash
# Two ten-step native VASP cell-relaxation pilots at 45.7 GPa from audited seeds.
# Outputs are diagnostic until ionic convergence, raw energies and phase identity pass.
#SBATCH -J varneb-gan-basin-pilot
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --array=0-1%2
#SBATCH --time=01:00:00

set -eo pipefail
source /public/home/iai806/Software/VASP/env.sh 6.3.2
set -u
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_NTASKS:-0}" == 32 ]] || exit 2
: "${WORK_ROOT:?set isolated WORK_ROOT}"
python=/public/home/iai806/.conda/envs/icu/bin/python
case_name=$("${python}" -c 'import json,sys; m=json.load(open(sys.argv[1])); assert m["status"]=="inputs_finalized_no_DFT" and m["purpose"]=="GaN_1000eV_signed_basin_native_VASP_10step_pilot_not_endpoint_certificate" and m["PSTRESS_kbar"]==457.0 and len(m["cases"])==2; print(m["cases"][int(sys.argv[2])]["name"])' "${WORK_ROOT}/manifest.json" "${SLURM_ARRAY_TASK_ID}")
case_dir="${WORK_ROOT}/cases/${case_name}"
[[ -s "${case_dir}/sha256.inputs.json" && ! -e "${case_dir}/OUTCAR" ]] || exit 3
"${python}" -c 'import hashlib,json,pathlib,sys; p=pathlib.Path(sys.argv[1]); d=json.loads((p/"sha256.inputs.json").read_text()); assert all(hashlib.sha256((p/k).read_bytes()).hexdigest()==v for k,v in d.items()); assert b"\r" not in (p/"INCAR").read_bytes(), "VASP INCAR must use LF only"' "${case_dir}"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
cd "${case_dir}"
srun --exclusive --exact --nodes=1 --ntasks=32 --ntasks-per-node=32 --cpus-per-task=1 \
  --mpi=pmi2 "${VASP_ROOT}/bin/vasp_std" > vasp.stdout 2>&1
grep -q 'General timing and accounting informations' OUTCAR
