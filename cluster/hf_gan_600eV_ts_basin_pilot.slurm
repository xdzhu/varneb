#!/bin/bash
# Two guarded 600-eV GaN native VASP pressure-relaxation pilots.
# A ten-step limit cannot certify either endpoint or the candidate saddle.
#SBATCH -J varneb-gan600-basin
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --array=0-1%2
#SBATCH --time=02:00:00

set -eo pipefail
source /public/home/iai806/Software/VASP/env.sh 6.3.2
set -u
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_NTASKS:-0}" == 32 ]] || exit 2
: "${WORK_ROOT:?set an isolated WORK_ROOT}"
python=/public/home/iai806/.conda/envs/icu/bin/python
case_name=$("${python}" - "${WORK_ROOT}/manifest.json" "${SLURM_ARRAY_TASK_ID}" <<'PY'
import json
import sys

record = json.load(open(sys.argv[1], encoding="utf-8"))
assert record["status"] == "inputs_finalized_no_DFT"
assert record["purpose"] == "GaN_45p7_600eV_signed_native_VASP_basin_10step_pilot"
assert record["encut_eV"] == 600 and record["PSTRESS_kbar"] == 457.0
assert record["ionic_settings"] == {"IBRION": 2, "ISIF": 3, "NSW": 10,
                                    "EDIFFG_eV_per_A": -0.02, "POTIM": 0.25}
assert len(record["cases"]) == 2
print(record["cases"][int(sys.argv[2])]["name"])
PY
)
case_dir="${WORK_ROOT}/cases/${case_name}"
[[ -s "${case_dir}/sha256.inputs.json" && ! -e "${case_dir}/OUTCAR" ]] || exit 3
"${python}" - "${case_dir}" <<'PY'
import hashlib
import json
import pathlib
import sys

directory = pathlib.Path(sys.argv[1])
expected = json.loads((directory / "sha256.inputs.json").read_text(encoding="utf-8"))
assert all(hashlib.sha256((directory / name).read_bytes()).hexdigest() == digest
           for name, digest in expected.items())
raw = (directory / "INCAR").read_bytes()
assert b"\r" not in raw and b"ENCUT = 600.000000\n" in raw
assert b"PSTRESS = 457.0\n" in raw and b"ISIF = 3\n" in raw
PY
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
cd "${case_dir}"
srun --exclusive --exact --nodes=1 --ntasks=32 --ntasks-per-node=32 \
  --cpus-per-task=1 --mpi=pmi2 "${VASP_ROOT}/bin/vasp_std" > vasp.stdout 2>&1
grep -q 'General timing and accounting informations' OUTCAR
