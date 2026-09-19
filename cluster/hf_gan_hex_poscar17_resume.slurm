#!/bin/bash
# Fresh FIRE from the last complete chain; no change to the physical contract.
#SBATCH -J gan-hex-poscar17
#SBATCH -p hfacnormal01
#SBATCH -N 3
#SBATCH --ntasks=291
#SBATCH --ntasks-per-node=97
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
set -eo pipefail
source /public/home/iai806/Software/VASP/env.sh 6.3.2
set -u
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_NTASKS}" == 291 && "${SLURM_JOB_NUM_NODES}" == 3 ]] || exit 2
: "${SOURCE_REPO:?}" "${OLD_TRAJECTORY:?}" "${CASE_ROOT:?}" "${NATIVE_PROBE:?}" "${WORKDIR:?}"
[[ ! -e "${WORKDIR}" ]] || { echo "workdir exists; refuse duplicate continuation" >&2; exit 3; }
python=/public/home/iai806/.conda/envs/icu/bin/python
export OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
export VASP_COMMAND="srun --exclusive --exact --nodes=1 --ntasks=32 --ntasks-per-node=32 --cpus-per-task=1 --distribution=block:block --mpi=pmi2 ${VASP_ROOT}/bin/vasp_std"
writer_sha=$(sha256sum "${SOURCE_REPO}/vcneb/vasp_contract.py" | cut -d' ' -f1)
vasp_sha=$(sha256sum "${VASP_ROOT}/bin/vasp_std" | cut -d' ' -f1)
export IMAGE_CACHE_NAMESPACE="hf-vasp6.3.2-${vasp_sha}-poscar17-${writer_sha}"
cd "${SOURCE_REPO}"
"${python}" examples/run_vcneb_vasp.py \
  --initial "${CASE_ROOT}/input/initial" --final "${CASE_ROOT}/input/final" \
  --resume --resume-trajectory "${OLD_TRAJECTORY}" \
  --workdir "${WORKDIR}" --n-images 29 --fmax 0.10 --steps 300 --k 0.20 --pressure-gpa 45.7 \
  --optimizer FIRE --image-workers 9 --no-climb \
  --initial-static-summary "${CASE_ROOT}/static_initial/vasp_static_summary.json" \
  --final-static-summary "${CASE_ROOT}/static_final/vasp_static_summary.json" \
  --image-cache-namespace "${IMAGE_CACHE_NAMESPACE}" \
  --native-lattice-probe "${NATIVE_PROBE}" \
  --candidate-step-retries 8 --candidate-step-retry-factor 0.5 \
  --minimum-distance 1.4 --maximum-deformation 0.50 --ncores 32 --vasp-bin "${VASP_ROOT}/bin/vasp_std"
"${python}" scripts/audit_vcneb_result.py "${WORKDIR}" \
  --max-min-distance 1.4 --max-deformation 0.50 --fmax-target 0.10
