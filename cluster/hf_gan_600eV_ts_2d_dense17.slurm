#!/bin/bash
# Only the 208 missing nodes of the nested 17x17 frozen joint-coordinate cut.
# This is the unchanged GaN VASP/PBE 600-eV, 45.7-GPa calculation contract.
#SBATCH -J varneb-gan600-17x17
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --array=0-207%4
#SBATCH --time=00:30:00
#SBATCH -o /public/home/iai806/abacus/agent-runs/gan_600eV_joint_17x17_20260930/slurm-%A_%a.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/gan_600eV_joint_17x17_20260930/slurm-%A_%a.err

set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_NTASKS:-0}" == 32 \
   && "${SLURM_JOB_PARTITION:-}" == hfacnormal01 ]] || exit 2
root=/public/home/iai806/abacus/agent-runs/gan_600eV_joint_17x17_20260930/new_cases
python=/public/home/iai806/.conda/envs/icu/bin/python
manifest_sha=ae913fdb7ecb656c7e4dbfa93e835028dd948e6d37acb919d94817cf8a35054b
echo "${manifest_sha}  ${root}/manifest.json" | sha256sum -c -
case_name=$("${python}" - "${root}/manifest.json" "${SLURM_ARRAY_TASK_ID}" <<'PY'
import json
import sys

record = json.load(open(sys.argv[1], encoding="utf-8"))
assert record["status"] == "inputs_finalized_no_DFT"
assert record["purpose"] == "GaN_45p7_600eV_local_joint_mode_17x17_frozen_refinement"
assert record["pressure_GPa"] == 45.7
assert record["n_prior_coordinates"] == 81 and len(record["cases"]) == 208
assert record["input_contract_sha256"]["INCAR"] == "83ba34d4b14ff4ea641bd74dec26e462ea1cc24b575db79a07138ccdfd552772"
print(record["cases"][int(sys.argv[2])]["name"])
PY
)
case_dir="${root}/cases/${case_name}"
[[ -s "${case_dir}/sha256.inputs.json" && ! -e "${case_dir}/OUTCAR" ]] || exit 3
source /public/home/iai806/Software/VASP/env.sh 6.3.2
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
