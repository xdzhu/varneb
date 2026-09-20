#!/bin/bash
# Generic Hefei VASP continuation for a preserved complete VCNEB chain.
# Resource count and all physical settings are supplied explicitly by sbatch.
# Endpoints are loaded from completed static summaries and are not recomputed.
#SBATCH -p hfacnormal01
#SBATCH --ntasks-per-node=97
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00

set -eo pipefail
source /public/home/iai806/Software/VASP/env.sh 6.3.2
set -u

required=(SOURCE_REPO INITIAL_DIR FINAL_DIR OLD_TRAJECTORY WORKDIR N_IMAGES
          IMAGE_WORKERS K_SPRING PRESSURE_GPA STATIC_INITIAL STATIC_FINAL
          CACHE_DIR CACHE_NAMESPACE MIN_DISTANCE MAX_DEFORMATION NATIVE_PROBE
          VASP_SYMPREC OPTIMIZER MAXSTEP)
for name in "${required[@]}"; do
  [[ -n "${!name:-}" ]] || { echo "missing required variable: ${name}" >&2; exit 2; }
done
[[ "${RUN_DFT:-0}" == 1 ]] || { echo "RUN_DFT=1 is required" >&2; exit 2; }
expected_tasks=$((SLURM_JOB_NUM_NODES * 97))
[[ "${SLURM_NTASKS}" == "${expected_tasks}" ]] || { echo "unexpected Slurm allocation" >&2; exit 2; }
[[ $((IMAGE_WORKERS * 32 + 1)) -le "${SLURM_NTASKS}" ]] || { echo "worker allocation exceeds task count" >&2; exit 2; }
[[ ! -e "${WORKDIR}" ]] || { echo "workdir exists; refusing duplicate continuation" >&2; exit 3; }
for path in "${SOURCE_REPO}" "${INITIAL_DIR}" "${FINAL_DIR}" "${OLD_TRAJECTORY}" \
            "${STATIC_INITIAL}" "${STATIC_FINAL}" "${NATIVE_PROBE}"; do
  [[ -e "${path}" ]] || { echo "missing input: ${path}" >&2; exit 2; }
done

python=/public/home/iai806/.conda/envs/icu/bin/python
export OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
export VASP_COMMAND="srun --exclusive --exact --nodes=1 --ntasks=32 --ntasks-per-node=32 --cpus-per-task=1 --distribution=block:block --mpi=pmi2 ${VASP_ROOT}/bin/vasp_std"
export IMAGE_CACHE_NAMESPACE="${CACHE_NAMESPACE}"
export VCNEB_GIT_REVISION="${SOURCE_REVISION:-unknown}"

optimizer_args=(--optimizer "${OPTIMIZER}" --maxstep "${MAXSTEP}")
candidate_args=(--native-lattice-probe "${NATIVE_PROBE}")
case "${OPTIMIZER}" in
  StagedFIRE)
    optimizer_args+=(--switch-fmax "${SWITCH_FMAX:-0.12}" --refine-maxstep "${REFINE_MAXSTEP:-0.02}")
    candidate_args+=(--candidate-step-retries "${CANDIDATE_STEP_RETRIES:-8}" \
                     --candidate-step-retry-factor "${CANDIDATE_STEP_RETRY_FACTOR:-0.5}")
    ;;
  FIRE|BlockFIRE|SplitFIRE)
    candidate_args+=(--candidate-step-retries "${CANDIDATE_STEP_RETRIES:-8}" \
                     --candidate-step-retry-factor "${CANDIDATE_STEP_RETRY_FACTOR:-0.5}")
    ;;
esac
if [[ "${OPTIMIZER}" == SplitFIRE ]]; then
  optimizer_args+=(--cell-maxstep "${CELL_MAXSTEP:-0.01}")
fi

mkdir -p "$(dirname "${WORKDIR}")"
cd "${SOURCE_REPO}"
# Optional retries repeat the identical calculator call; they never alter the contract.
"${python}" examples/run_vcneb_vasp.py \
  --initial "${INITIAL_DIR}" --final "${FINAL_DIR}" \
  --resume --resume-trajectory "${OLD_TRAJECTORY}" \
  --workdir "${WORKDIR}" --n-images "${N_IMAGES}" --fmax 0.10 \
  --steps "${STEPS:-300}" --k "${K_SPRING}" --pressure-gpa "${PRESSURE_GPA}" \
  "${optimizer_args[@]}" --image-workers "${IMAGE_WORKERS}" \
  --image-retries "${IMAGE_RETRIES:-0}" --no-climb \
  --initial-static-summary "${STATIC_INITIAL}" --final-static-summary "${STATIC_FINAL}" \
  --image-cache-dir "${CACHE_DIR}" --image-cache-namespace "${IMAGE_CACHE_NAMESPACE}" \
  --vasp-isym -1 --vasp-symprec "${VASP_SYMPREC}" "${candidate_args[@]}" \
  --minimum-distance "${MIN_DISTANCE}" --maximum-deformation "${MAX_DEFORMATION}" \
  --ncores 32 --vasp-bin "${VASP_ROOT}/bin/vasp_std"

"${python}" scripts/audit_vcneb_result.py "${WORKDIR}" \
  --max-min-distance "${MIN_DISTANCE}" --max-deformation "${MAX_DEFORMATION}" --fmax-target 0.10
