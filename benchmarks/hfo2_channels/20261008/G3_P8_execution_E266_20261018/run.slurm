#!/bin/bash
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=04:00:00
set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_JOB_PARTITION:-}" == hfacnormal01 \
 && "${SLURM_JOB_NUM_NODES:-0}" == 1 && "${SLURM_NTASKS:-0}" == 1 \
 && "${SLURM_CPUS_PER_TASK:-0}" == 32 && "${SHARD:-x}" =~ ^[01]$ \
 && "${SLURM_JOB_NAME:-}" == "hfo2-G3-E266-P${SHARD}" ]] || exit 2
bundle=${BUNDLE:?immutable bundle required}
/public/home/iai806/.conda/envs/icu/bin/python -c '
import hashlib,json,pathlib,sys
p=pathlib.Path(sys.argv[1]);m=p/"bundle_files.json"
if hashlib.sha256(m.read_bytes()).hexdigest()!=sys.argv[2]:raise ValueError("bundle manifest drift")
files=json.loads(m.read_bytes())
expected={"worker.py","run.slurm","config.json","forecasts.json","starts.json","protocol.md"}
if set(files)!=expected or any(pathlib.PurePosixPath(n).name!=n or hashlib.sha256((p/n).read_bytes()).hexdigest()!=h for n,h in files.items()):raise ValueError("bundle drift")
' "${bundle}" "${BUNDLE_MANIFEST_SHA256:?exact bundle hash required}"
module purge
source /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS
export I_MPI_HYDRA_BOOTSTRAP=fork
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH=/public/home/iai806/abacus/agent-runs/20261015-varneb-G3-faces-E161-r1/source
cd "${PYTHONPATH}"
/public/home/iai806/.conda/envs/icu/bin/python "${bundle}/worker.py" run \
 --bundle "${bundle}" --config "${bundle}/config.json" \
 --config-sha256 "${CONFIG_SHA256:?exact config required}" --shard "${SHARD}"
