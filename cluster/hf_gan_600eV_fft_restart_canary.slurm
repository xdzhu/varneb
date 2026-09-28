#!/bin/bash
# Isolated one-step VASP test of the original GaN pilot FFT grid.
#SBATCH -J varneb-gan600-fft
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00

set -eo pipefail
source /public/home/iai806/Software/VASP/env.sh 6.3.2
set -u
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_NTASKS:-0}" == 32 ]] || exit 2
: "${WORK_ROOT:?set an isolated WORK_ROOT}"
python=/public/home/iai806/.conda/envs/icu/bin/python
"${python}" - "${WORK_ROOT}" <<'PY'
import hashlib
import json
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
record = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
assert record["status"] == "inputs_finalized_no_DFT"
assert record["purpose"] == "GaN_45p7_600eV_negative_basin_original_FFT_grid_one_step_canary"
assert record["encut_eV"] == 600 and record["PSTRESS_kbar"] == 457.0
directory = root / "case"
assert not (directory / "OUTCAR").exists()
assert all(hashlib.sha256((directory / name).read_bytes()).hexdigest() == digest
           for name, digest in record["input_sha256"].items())
raw = (directory / "INCAR").read_bytes()
assert b"\r" not in raw and b"ENCUT = 600.000000\n" in raw
assert b"PSTRESS = 457.0\n" in raw and b"NGZ = 40\n" in raw
PY
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
cd "${WORK_ROOT}/case"
srun --exclusive --exact --nodes=1 --ntasks=32 --ntasks-per-node=32 \
  --cpus-per-task=1 --mpi=pmi2 "${VASP_ROOT}/bin/vasp_std" > vasp.stdout 2>&1
grep -q 'General timing and accounting informations' OUTCAR
