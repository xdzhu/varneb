#!/bin/bash
#SBATCH -J varneb-gan-path-check
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=03:00:00
set -eo pipefail
source /public/home/iai806/Software/VASP/env.sh 6.3.2
set -u
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_NTASKS}" == 32 ]] || exit 2
repo=/public/home/iai806/abacus/agent-runs/20260918-varneb-gan-suite
parent=/public/home/iai806/abacus/agent-runs/20260918-varneb-vasp-hf/gan_b4_b1_vasp_pbe_paw_qian
python=/public/home/iai806/.conda/envs/icu/bin/python
export OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
export VASP_COMMAND="srun --exclusive --exact --nodes=1 --ntasks=32 --ntasks-per-node=32 --cpus-per-task=1 --mpi=pmi2 ${VASP_ROOT}/bin/vasp_std"
cd "${repo}"
case "${CASE_NAME}" in
  b4_b1_hexagonal) pressure=45.7 ;;
  b3_b1_diagonal)
    pressure=45.0
    "${python}" scripts/prepare_gan_qian_suite.py --stage b3 --parent-case "${parent}" --suite-root cases
    ;;
  *) exit 2 ;;
esac
case_root="${repo}/cases/${CASE_NAME}"
"${python}" examples/run_vcneb_vasp.py \
  --initial "${case_root}/input/initial" --final "${case_root}/input/final" \
  --initial-trajectory "${case_root}/initial.traj" --workdir "${case_root}/preflight" \
  --n-images 29 --pressure-gpa "${pressure}" --no-climb --validate-only \
  --minimum-distance 1.4 --maximum-deformation 0.50 --ncores 32 --vasp-bin "${VASP_ROOT}/bin/vasp_std"
"${python}" scripts/probe_vasp_input_contract.py \
  --source-dir "${case_root}/input/initial" --trajectory "${case_root}/initial.traj" \
  --n-images 29 --workdir "${case_root}/initialization_probe" --symprec 1e-4 \
  --ncores 32 --vasp-bin "${VASP_ROOT}/bin/vasp_std" --timeout 180 --run
for endpoint in initial final; do
  "${python}" examples/run_vcneb_vasp.py \
    --initial "${case_root}/input/initial" --final "${case_root}/input/final" \
    --initial-trajectory "${case_root}/initial.traj" --workdir "${case_root}/static_${endpoint}" \
    --n-images 29 --pressure-gpa "${pressure}" --no-climb --static-only --static-endpoint "${endpoint}" \
    --minimum-distance 1.4 --maximum-deformation 0.50 --ncores 32 --vasp-bin "${VASP_ROOT}/bin/vasp_std"
done
