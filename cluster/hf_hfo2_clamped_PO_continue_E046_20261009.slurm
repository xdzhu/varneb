#!/bin/bash
# One reviewed PO+ geometry continuation; no automatic matrix or restart.
#SBATCH -J hfo2-PO-E046
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=02:00:00
set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_JOB_NUM_NODES:-0}" == 1 \
   && "${SLURM_NTASKS:-0}" == 1 && "${SLURM_CPUS_PER_TASK:-0}" == 32 ]] || exit 2
source_root=${SOURCE_ROOT:?immutable tested source required}
endpoint_root=${ENDPOINT_ROOT:?fresh endpoint namespace required}
seed_root=${SEED_ROOT:?fresh HF full-byte audited continuation seed required}
parent=/public/home/iai806/abacus/agent-runs/20261009-varneb-clamped-canary-E045-r1/endpoint
[[ ! -e "${endpoint_root}" && -f "${seed_root}/preparation.json" ]] || exit 2
module purge
source /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS
export I_MPI_HYDRA_BOOTSTRAP=fork
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="${source_root}${PYTHONPATH:+:${PYTHONPATH}}"
cd "${source_root}"
python=/public/home/iai806/.conda/envs/icu/bin/python
"${python}" -c 'import json,sys; p=json.load(open(sys.argv[1])); assert p["new_DFT_calls"]==0 and p["parent_review_full_physical_bytes_checked_here"] is True and p["planned_max_BFGS_steps"]==p["planned_max_new_SCF_calls"]==20 and p["identical_initial_SCF_reused"] is True' "${seed_root}/preparation.json"
set +e
"${python}" -m scripts.relax_clamped_ase_endpoint --seed-manifest "${seed_root}/endpoint_seed.json" \
  --workdir "${endpoint_root}" --factory examples.hfo2_fixed_input_factory:make_endpoint_cached_factory \
  --parameters "${seed_root}/factory_parameters.json" \
  --command "mpirun -np 32 /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus" \
  --fmax 0.03 --stress-kbar 2 --maxstep 0.02 --steps 20
relax_exit=$?
set -e
[[ ${relax_exit} == 0 || ${relax_exit} == 1 ]] || exit "${relax_exit}"
"${python}" -m scripts.audit_hfo2_clamped_continuation --root "${endpoint_root}" \
  --seed-manifest "${seed_root}/endpoint_seed.json" --parent-endpoint "${parent}" \
  --output "${endpoint_root}/continuation_audit.json"
