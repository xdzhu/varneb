#!/bin/bash
# One manually admitted bounded treatment ablation, not a monitor.
#SBATCH -J hfo2-G3-E175-atoms
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=05:00:00
set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_JOB_PARTITION:-}" == hfacnormal01 \
 && "${SLURM_JOB_NUM_NODES:-0}" == 1 && "${SLURM_NTASKS:-0}" == 1 \
 && "${SLURM_CPUS_PER_TASK:-0}" == 32 ]] || exit 2
bundle=${BUNDLE:?immutable external bundle required}
bundle_sha=${BUNDLE_MANIFEST_SHA256:?exact bundle manifest required}
config_sha=${CONFIG_SHA256:?exact atomic-release configuration required}
/public/home/iai806/.conda/envs/icu/bin/python -c '
import hashlib,json,pathlib,sys
p=pathlib.Path(sys.argv[1]);m=p/"bundle_files.json"
if hashlib.sha256(m.read_bytes()).hexdigest()!=sys.argv[2]: raise ValueError("bundle manifest drift")
files=json.loads(m.read_bytes())
if len(files)!=5 or any(pathlib.PurePosixPath(n).name!=n or hashlib.sha256((p/n).read_bytes()).hexdigest()!=h for n,h in files.items()): raise ValueError("bundle files drift")
' "${bundle}" "${bundle_sha}"
module purge
source /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS
export I_MPI_HYDRA_BOOTSTRAP=fork
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH=/public/home/iai806/abacus/agent-runs/20261015-varneb-G3-faces-E161-r1/source
cd "${PYTHONPATH}"
/public/home/iai806/.conda/envs/icu/bin/python "${bundle}/worker.py" run \
 --bundle "${bundle}" --config "${bundle}/config.json" --config-sha256 "${config_sha}" \
 --workdir "${bundle}/../atomic_release"
