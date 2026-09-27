#!/bin/bash
# Bounded GaN central-path atomic-transverse tube; original VASP 600 eV.
#SBATCH -J varneb-gan600-atomic-central
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --array=0-15%4
#SBATCH --time=00:30:00

set -eo pipefail
source /public/home/iai806/Software/VASP/env.sh 6.3.2
set -u
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_NTASKS:-0}" == 32 ]] || exit 2
: "${WORK_ROOT:?set the finalized GaN central atomic-tube directory}"
python=/public/home/iai806/.conda/envs/icu/bin/python
case_name=$("${python}" - "${WORK_ROOT}/manifest.json" "${SLURM_ARRAY_TASK_ID}" <<'PY'
import json
import sys

record = json.load(open(sys.argv[1], encoding="utf-8"))
assert record["status"] == "inputs_finalized_no_DFT"
assert record["purpose"] == "GaN_45p7_600eV_atomic_only_central_path_tube_pilot"
assert record["anchors"] == [5, 8, 11, 14, 17, 18, 20, 22]
assert record["q_atom_A"] == [-0.05, 0.0, 0.05]
assert record["predeclared_LOO_gate_meV_per_GaN"] == 1.0
assert record["input_contract_sha256"]["INCAR"] == "83ba34d4b14ff4ea641bd74dec26e462ea1cc24b575db79a07138ccdfd552772"
assert len(record["cases"]) == 16
print(record["cases"][int(sys.argv[2])]["name"])
PY
)
case_dir="${WORK_ROOT}/cases/${case_name}"
[[ -s "${case_dir}/sha256.inputs.json" && ! -e "${case_dir}/OUTCAR" ]] || exit 3
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
cd "${case_dir}"
"${python}" - <<'PY'
import hashlib
import json
from pathlib import Path

expected = json.loads(Path("sha256.inputs.json").read_text(encoding="utf-8"))
assert all(hashlib.sha256(Path(name).read_bytes()).hexdigest() == digest
           for name, digest in expected.items())
assert expected["INCAR"] == "83ba34d4b14ff4ea641bd74dec26e462ea1cc24b575db79a07138ccdfd552772"
PY
printf 'job=%s task=%s node=%s case=%s\n' \
  "${SLURM_JOB_ID}" "${SLURM_ARRAY_TASK_ID}" "${SLURM_JOB_NODELIST}" "${case_name}"
srun --exclusive --exact --nodes=1 --ntasks=32 --ntasks-per-node=32 --cpus-per-task=1 \
  --mpi=pmi2 "${VASP_ROOT}/bin/vasp_std" > vasp.stdout 2>&1
grep -q 'General timing and accounting informations' OUTCAR
