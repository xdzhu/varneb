#!/bin/bash
# Bounded 12-case test of three normal-strain derivatives at unchanged 600 eV.
# Submit explicitly with WORK_ROOT and RUN_DFT=1; never a CI job.
#SBATCH -J varneb-gan600-dHdF
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --array=0-11%4
#SBATCH --time=00:30:00

set -eo pipefail
source /public/home/iai806/Software/VASP/env.sh 6.3.2
set -u
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_JOB_PARTITION:-}" == hfacnormal01 \
   && "${SLURM_NTASKS:-0}" == 32 ]] || exit 2
: "${WORK_ROOT:?set the isolated prepared work root}"
python=/public/home/iai806/.conda/envs/icu/bin/python
case_name=$("${python}" - "${WORK_ROOT}/manifest.json" "${SLURM_ARRAY_TASK_ID}" <<'PY'
import json
import sys
from pathlib import Path

report = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
assert report["status"] == "inputs_finalized_no_DFT"
assert report["purpose"] == "GaN_600eV_normal_strain_derivative_step_test"
assert report["pressure_GPa"] == 45.7
assert report["n_cases"] == 12
assert report["electronic_input_sha256"]["INCAR"] == "83ba34d4b14ff4ea641bd74dec26e462ea1cc24b575db79a07138ccdfd552772"
assert report["electronic_input_sha256"]["POTCAR"] == "f94781ce6cf9b9454c093353320c5e80a2a0aceea6fb8353b33b01ac37b95168"
case = report["cases"][int(sys.argv[2])]
assert case["step_A"] in (0.01, 0.005) and case["axis"] in (12, 13, 14)
print(case["name"])
PY
)
case_dir="${WORK_ROOT}/cases/${case_name}"
[[ -s "${case_dir}/sha256.inputs.json" && ! -e "${case_dir}/OUTCAR" ]] || exit 3
"${python}" - "${case_dir}" <<'PY'
import hashlib
import json
import sys
from pathlib import Path

directory = Path(sys.argv[1])
expected = json.loads((directory / "sha256.inputs.json").read_text())
assert set(expected) == {"INCAR", "KPOINTS", "POSCAR", "POTCAR"}
assert all(hashlib.sha256((directory / key).read_bytes()).hexdigest() == value
           for key, value in expected.items())
PY
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
cd "${case_dir}"
srun --exclusive --exact --nodes=1 --ntasks=32 --ntasks-per-node=32 \
  --cpus-per-task=1 --mpi=pmi2 "${VASP_ROOT}/bin/vasp_std" > vasp.stdout 2>&1
grep -q 'General timing and accounting informations' OUTCAR
