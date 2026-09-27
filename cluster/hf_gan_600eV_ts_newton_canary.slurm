#!/bin/bash
# One GaN 45.7-GPa image-15 Newton candidate: original 600-eV VASP contract.
# This is a static test, not a claimed saddle or a new VCNEB production chain.
#SBATCH -J varneb-gan600-ts1
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00

set -eo pipefail
source /public/home/iai806/Software/VASP/env.sh 6.3.2
set -u
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_NTASKS:-0}" == 32 ]] || exit 2
: "${WORK_ROOT:?set the staged GaN 600-eV canary directory}"
[[ -s "${WORK_ROOT}/manifest.json" && ! -e "${WORK_ROOT}/OUTCAR" ]] || exit 3
python=/public/home/iai806/.conda/envs/icu/bin/python
"${python}" - "${WORK_ROOT}" <<'PY'
import hashlib
import json
from pathlib import Path
import sys

root = Path(sys.argv[1])
record = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
assert record["status"] == "geometry_preflight_passed_no_DFT"
assert record["purpose"] == "GaN_45p7_600eV_image15_one_step_joint_Newton_static_canary"
assert record["electronic_contract"]["ENCUT_eV"] == 600
expected = dict(record["electronic_contract"]["original_input_sha256"])
expected["POSCAR"] = record["POSCAR_sha256"]
for filename, digest in expected.items():
    actual = hashlib.sha256((root / filename).read_bytes()).hexdigest()
    assert actual == digest, (filename, actual, digest)
assert not (root / "OUTCAR").exists()
PY
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
cd "${WORK_ROOT}"
printf 'job=%s node=%s work=%s\n' "${SLURM_JOB_ID}" "${SLURM_JOB_NODELIST}" "${WORK_ROOT}"
srun --exclusive --exact --nodes=1 --ntasks=32 --ntasks-per-node=32 --cpus-per-task=1 \
  --mpi=pmi2 "${VASP_ROOT}/bin/vasp_std" > vasp.stdout 2>&1
grep -q 'General timing and accounting informations' OUTCAR
