#!/bin/bash
# Exactly one reviewed zero-strain PO+ canary, not a production matrix.
#SBATCH -J hfo2-clamp-E045
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=00:30:00
set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_JOB_NUM_NODES:-0}" == 1 \
   && "${SLURM_NTASKS:-0}" == 1 && "${SLURM_CPUS_PER_TASK:-0}" == 32 ]] || exit 2
source_root=${SOURCE_ROOT:?immutable tested source required}
endpoint_root=${ENDPOINT_ROOT:?fresh endpoint namespace required}
gate=${G1_REVIEW:?fresh HF raw G1 candidate review required}
[[ ! -e "${endpoint_root}" && -f "${gate}" ]] || exit 2
module purge
source /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS
export I_MPI_HYDRA_BOOTSTRAP=fork
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="${source_root}${PYTHONPATH:+:${PYTHONPATH}}"
cd "${source_root}"
python=/public/home/iai806/.conda/envs/icu/bin/python
"${python}" -c 'import json,sys; g=json.load(open(sys.argv[1])); assert g["candidate_G1_gate_passed"] is True and g["raw_on_HF"] is True and g["fresh_raw_terminal_records"]==40 and g["new_DFT_calls"]==0' "${gate}"
seed="${source_root}/benchmarks/hfo2_channels/20261008/clamped_endpoint_seeds/strain_0000/PO_plus/endpoint_seed.json"
parameters="${source_root}/benchmarks/hfo2_channels/20261008/M_endpoint_factory_parameters.json"
set +e
"${python}" -m scripts.relax_clamped_ase_endpoint --seed-manifest "${seed}" --workdir "${endpoint_root}" \
  --factory examples.hfo2_fixed_input_factory:make_factory --parameters "${parameters}" \
  --command "mpirun -np 32 /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus" \
  --fmax 0.03 --stress-kbar 2 --maxstep 0.02 --steps 4
relax_exit=$?
set -e
[[ ${relax_exit} == 0 || ${relax_exit} == 1 ]] || exit "${relax_exit}"
# Exit1 is accepted only if a complete, raw-replayed step-limit summary exists.
# Unexpected failures without that evidence still make this job fail.
"${python}" -m scripts.audit_hfo2_clamped_canary --root "${endpoint_root}" \
  --seed-manifest "${seed}" --output "${endpoint_root}/canary_audit.json"
