#!/bin/bash
# 4-atom GaN highest-image local joint curvature; no phonon supercell.
#SBATCH -J varneb-gan-jointH
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --array=0-35%4
#SBATCH --time=01:00:00

set -eo pipefail
source /public/home/iai806/Software/VASP/env.sh 6.3.2
set -u
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_NTASKS:-0}" == 32 ]] || exit 2
: "${WORK_ROOT:?set immutable staged WORK_ROOT}"
[[ -s "${WORK_ROOT}/manifest.json" ]] || exit 2
python=/public/home/iai806/.conda/envs/icu/bin/python
case_name=$("${python}" -c 'import json,sys; d=json.load(open(sys.argv[1])); assert d["status"]=="inputs_finalized_no_DFT"; print(d["cases"][int(sys.argv[2])]["name"])' \
  "${WORK_ROOT}/manifest.json" "${SLURM_ARRAY_TASK_ID}")
case_dir="${WORK_ROOT}/cases/${case_name}"
[[ -s "${case_dir}/sha256.inputs.json" && ! -e "${case_dir}/OUTCAR" ]] || exit 3
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
cd "${case_dir}"
"${python}" -c 'import hashlib,json,pathlib; d=json.load(open("sha256.inputs.json")); assert all(hashlib.sha256(pathlib.Path(k).read_bytes()).hexdigest()==v for k,v in d.items())'
printf 'job=%s task=%s node=%s case=%s\n' \
  "${SLURM_JOB_ID}" "${SLURM_ARRAY_TASK_ID}" "${SLURM_JOB_NODELIST}" "${case_name}"
srun --exclusive --exact --nodes=1 --ntasks=32 --ntasks-per-node=32 --cpus-per-task=1 \
  --mpi=pmi2 "${VASP_ROOT}/bin/vasp_std" > vasp.stdout 2>&1
grep -q 'General timing and accounting informations' OUTCAR
