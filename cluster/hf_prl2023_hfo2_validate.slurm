#!/bin/bash
#SBATCH -J hfo2-prl-check
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
set -eo pipefail
source /public/home/iai806/Software/VASP/env.sh 6.3.2
set -u

[[ "${RUN_DFT:-0}" == 1 && "${SLURM_NTASKS}" == 32 ]] || exit 2
case "${CASE_NAME:-}" in t_to_po|po_to_m) ;; *) exit 2 ;; esac

repo=/public/home/iai806/abacus/agent-runs/20260919-prl2023-hfo2-fig2a
python=/public/home/iai806/.conda/envs/icu/bin/python
case_root="${repo}/cases/${CASE_NAME}"
initial="${case_root}/endpoint_relax/initial"
final="${case_root}/endpoint_relax/final"
path="${case_root}/initial_corrected.traj"
export OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
export VASP_COMMAND="srun --exclusive --exact --nodes=1 --ntasks=32 --ntasks-per-node=32 --cpus-per-task=1 --mpi=pmi2 ${VASP_ROOT}/bin/vasp_std"
cd "${repo}"

"${python}" scripts/prepare_prl2023_hfo2_fig2a.py correct \
  --case-root "${case_root}" --initial "${initial}/CONTCAR" --final "${final}/CONTCAR" --output "${path}"
"${python}" examples/run_vcneb_vasp.py \
  --initial "${initial}" --final "${final}" --initial-trajectory "${path}" \
  --workdir "${case_root}/preflight" --n-images 20 --fmax 0.10 --steps 300 --k 0.10 \
  --pressure-gpa 0 --optimizer FIRE --no-climb --validate-only \
  --minimum-distance 1.40 --maximum-deformation 0.50 --ncores 32 --vasp-bin "${VASP_ROOT}/bin/vasp_std"
"${python}" - "${path}" "${case_root}/midpoint.vasp" <<'PY'
from ase.io import read, write
import sys
images = read(sys.argv[1], index=":")
write(sys.argv[2], images[len(images) // 2], format="vasp", direct=True, vasp5=True)
PY
"${python}" examples/run_vasp_single_image_smoke.py \
  --structure "${case_root}/midpoint.vasp" --template-dir "${initial}" \
  --workdir "${case_root}/midpoint_static" --ncores 32 --vasp-bin "${VASP_ROOT}/bin/vasp_std"
