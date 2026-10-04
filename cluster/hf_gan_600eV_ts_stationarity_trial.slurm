#!/bin/bash
# Nine independent, exact-contract static probes around one energy-informed
# GaN center. This is a stationarity test, not a new VCNEB production path.
#SBATCH -J gan600-ts-cert
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --array=0-8%3
#SBATCH --time=04:00:00
#SBATCH -o slurm-%x-%A_%a.out
#SBATCH -e slurm-%x-%A_%a.err
set -euo pipefail
source /public/home/iai806/Software/VASP/env.sh 6.3.2
[[ "${RUN_DFT:-0}" == 1 ]] || exit 2
: "${CASE_ROOT:?set the prepared stationarity-trial case directory}"
[[ "${SLURM_NTASKS:-0}" == 32 ]] || exit 2
python=/public/home/iai806/.conda/envs/icu/bin/python
case_dir="$("${python}" -c '
import hashlib, json, sys
from pathlib import Path
root, index = Path(sys.argv[1]), int(sys.argv[2])
m = json.loads((root / "manifest.json").read_text())
if m["status"] != "predeclared_stationarity_trial_inputs_ready_no_DFT" or m["n_cases"] != 9:
    raise SystemExit("invalid trial manifest")
case = m["cases"][index]
d = root / "cases" / case["name"]
for filename, digest in case["input_sha256"].items():
    if hashlib.sha256((d / filename).read_bytes()).hexdigest() != digest:
        raise SystemExit("input hash mismatch: " + case["name"] + "/" + filename)
print(d)
' "${CASE_ROOT}" "${SLURM_ARRAY_TASK_ID}")"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
cd "${case_dir}"
[[ ! -e OUTCAR ]] || exit 3
srun --mpi=pmi2 --ntasks=32 "${VASP_ROOT}/bin/vasp_std" > vasp.stdout
grep -q 'aborting loop because EDIFF is reached' OUTCAR
grep -q 'General timing and accounting informations for this job' OUTCAR
