#!/bin/bash
# Same-setting 1000 eV static check of an audited GaN basin candidate.
#SBATCH -J varneb-gan-basin-static
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=00:45:00

set -eo pipefail
source /public/home/iai806/Software/VASP/env.sh 6.3.2
set -u
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_NTASKS:-0}" == 32 ]] || exit 2
: "${WORK_ROOT:?set isolated WORK_ROOT}"
python=/public/home/iai806/.conda/envs/icu/bin/python
"${python}" -c 'import json,sys; m=json.load(open(sys.argv[1])); assert m["status"]=="inputs_finalized_no_DFT" and m["purpose"]=="GaN_1000eV_audited_basin_static_E_plus_PV_not_TS_certificate" and m["calculator"]["ENCUT_eV"]==1000' "${WORK_ROOT}/manifest.json"
case_dir="${WORK_ROOT}/case"
[[ -s "${case_dir}/sha256.inputs.json" && ! -e "${case_dir}/OUTCAR" ]] || exit 3
"${python}" -c 'import hashlib,json,pathlib,sys; p=pathlib.Path(sys.argv[1]); d=json.loads((p/"sha256.inputs.json").read_text()); assert all(hashlib.sha256((p/k).read_bytes()).hexdigest()==v for k,v in d.items()); assert b"\r" not in (p/"INCAR").read_bytes()' "${case_dir}"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
cd "${case_dir}"
srun --exclusive --exact --nodes=1 --ntasks=32 --ntasks-per-node=32 --cpus-per-task=1 \
  --mpi=pmi2 "${VASP_ROOT}/bin/vasp_std" > vasp.stdout 2>&1
grep -q 'General timing and accounting informations' OUTCAR
