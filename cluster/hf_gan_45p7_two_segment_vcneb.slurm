#!/bin/bash
# One of the two fixed-center GaN segments. Submit left and right separately.
# Each hfacnormal01 node hosts three 32-MPI workers plus manager headroom.
#SBATCH -J gan600-split
#SBATCH -p hfacnormal01
#SBATCH -N 2
#SBATCH --ntasks=194
#SBATCH --ntasks-per-node=97
#SBATCH --cpus-per-task=1
#SBATCH --time=96:00:00
#SBATCH -o slurm-%x-%j.out
#SBATCH -e slurm-%x-%j.err
set -euo pipefail
source /public/home/iai806/Software/VASP/env.sh 6.3.2
[[ "${RUN_DFT:-0}" == 1 ]] || exit 2
[[ "${SLURM_NTASKS:-0}" == 194 && "${SLURM_JOB_NUM_NODES:-0}" == 2 ]] || exit 2
: "${SPLIT_ROOT:?set the prepared two-segment directory}"
: "${CASE_NAME:?set the prepared case directory name}"
: "${SEGMENT:?set left or right}"
[[ "${SEGMENT}" == left || "${SEGMENT}" == right ]] || exit 2
[[ -s "${SPLIT_ROOT}/${CASE_NAME}/manifest.json" ]] || exit 3
[[ ! -e "${SPLIT_ROOT}/${CASE_NAME}/${SEGMENT}/run/vcneb.traj" ]] || exit 3

python=/public/home/iai806/.conda/envs/icu/bin/python
repo="${SPLIT_ROOT}/source"
case_dir="${SPLIT_ROOT}/${CASE_NAME}/${SEGMENT}"
if [[ "${SEGMENT}" == left ]]; then n_images=16; else n_images=14; fi
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONUNBUFFERED=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
export VASP_COMMAND="srun --exclusive --exact --nodes=1 --ntasks=32 --ntasks-per-node=32 --cpus-per-task=1 --distribution=block:block --mpi=pmi2 ${VASP_ROOT}/bin/vasp_std"
export VCNEB_GIT_REVISION="${VCNEB_GIT_REVISION:?pass the archived source commit}"
export IMAGE_CACHE_NAMESPACE="hf-vasp6.3.2-$(sha256sum "${VASP_ROOT}/bin/vasp_std" | cut -d' ' -f1)"
cd "${repo}"
"${python}" examples/run_vcneb_vasp.py \
  --initial "${case_dir}/initial" --final "${case_dir}/final" \
  --workdir "${case_dir}/run" --n-images "${n_images}" \
  --fmax 0.10 --steps 300 --k 0.20 --pressure-gpa 45.7 \
  --optimizer FIRE --image-workers 6 --no-climb \
  --initial-static-summary "${case_dir}/static_initial.json" \
  --final-static-summary "${case_dir}/static_final.json" \
  --resume --resume-trajectory "${case_dir}/seed.traj" \
  --image-cache-namespace "${IMAGE_CACHE_NAMESPACE}" \
  --mic --cell-interpolation linear --mapping identity \
  --minimum-distance 1.4 --maximum-deformation 0.50 \
  --vasp-isym -1 --vasp-symprec 1e-4 --expected-encut-ev 600 \
  --ncores 32 --vasp-bin "${VASP_ROOT}/bin/vasp_std"
"${python}" scripts/audit_vcneb_result.py "${case_dir}/run" \
  --max-min-distance 1.4 --max-deformation 0.50 --fmax-target 0.10
