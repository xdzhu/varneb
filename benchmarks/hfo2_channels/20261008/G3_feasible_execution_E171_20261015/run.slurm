#!/bin/bash
# Manual admission only, after completed/audited E166 controls. No timer.
#SBATCH -J hfo2-G3-E171-feasible
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00
set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_JOB_PARTITION:-}" == hfacnormal01 \
 && "${SLURM_JOB_NUM_NODES:-0}" == 1 && "${SLURM_NTASKS:-0}" == 1 \
 && "${SLURM_CPUS_PER_TASK:-0}" == 32 ]] || exit 2
bundle=${BUNDLE:?immutable external bundle required}
bundle_manifest_sha=${BUNDLE_MANIFEST_SHA256:?exact bundle manifest required}
config_sha=${EXECUTION_CONFIG_SHA256:?exact execution configuration required}
old=/public/home/iai806/abacus/agent-runs/20261015-varneb-G3-faces-E161-r1/source
/public/home/iai806/.conda/envs/icu/bin/python -c '
import hashlib,json,pathlib,sys
p=pathlib.Path(sys.argv[1]);m=p/"bundle_files.json"
if hashlib.sha256(m.read_bytes()).hexdigest()!=sys.argv[2]: raise ValueError("bundle manifest drift")
files=json.loads(m.read_bytes())
if not files or any(pathlib.PurePosixPath(n).name!=n or hashlib.sha256((p/n).read_bytes()).hexdigest()!=h for n,h in files.items()): raise ValueError("bundle file drift")
' "${bundle}" "${bundle_manifest_sha}"
module purge
source /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS
export I_MPI_HYDRA_BOOTSTRAP=fork
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="${old}"
cd "${old}"
/public/home/iai806/.conda/envs/icu/bin/python "${bundle}/worker.py" run \
 --bundle "${bundle}" --config "${bundle}/execution.json" --config-sha256 "${config_sha}" \
 --workdir "${bundle}/../arm"
