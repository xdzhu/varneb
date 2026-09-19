#!/bin/bash
# Resume the GaN B4->B1 serial VCNEB after a VASP Bravais-lattice refusal.
# The failed directory is read-only evidence; this creates a new recovery
# directory with a fixed symmetry policy. Probe the entire saved chain before
# endpoint statics or production VCNEB; never tune parameters per image.

set -euo pipefail

repo=${REPO:-/home/zhuxd/abacus/agent-runs/20260916-varneb-v}
python=${PYTHON:-/home/zhuxd/Software/anaconda3/envs/icu/bin/python}
vasp_bin=${VASP_BIN:-/home/zhuxd/Software/src/vasp/6.3.2/bin/vasp_std}
case_root=${CASE_ROOT:-${repo}/gan_b4_b1_vasp_pbe_paw_qian}
old_workdir=${OLD_WORKDIR:-${case_root}/vcneb_b4_to_b1_tetragonal_n29_cu17_serial}
recovery_root=${RECOVERY_ROOT:-${case_root}/recovery_isym_minus1}
pressure_gpa=${PRESSURE_GPA:-45.7}
vasp_symprec=${VASP_SYMPREC:-1e-4}
steps=${STEPS:-300}
fmax=${FMAX:-0.10}
spring=${SPRING:-0.20}

[[ "${RUN_DFT:-0}" == 1 ]] || { echo "Set RUN_DFT=1 to run VASP on cu17" >&2; exit 2; }
[[ "$(hostname -s)" == cu17 ]] || { echo "This runner is restricted to cu17" >&2; exit 2; }
# GNU nproc honors OpenMP thread limits; query CPU affinity without confusing
# OMP_NUM_THREADS=1 (one thread per MPI rank) with a one-core allocation.
visible_cores=$(env -u OMP_NUM_THREADS -u OMP_THREAD_LIMIT nproc)
[[ "${visible_cores}" == 40 ]] || { echo "cu17 serial workflow requires exactly 40 visible CPU cores" >&2; exit 2; }
: "${VCNEB_GIT_REVISION:?Set the exact source commit in VCNEB_GIT_REVISION}"
[[ -f "${old_workdir}/vcneb.traj" ]] || { echo "missing recovery trajectory: ${old_workdir}/vcneb.traj" >&2; exit 2; }
[[ ! -e "${recovery_root}" ]] || { echo "recovery root already exists: ${recovery_root}" >&2; exit 2; }

export OMP_NUM_THREADS=1
export VASP_COMMAND="mpirun -np 40 ${vasp_bin}"
cd "${repo}"
mkdir -p "${recovery_root}"

"${python}" scripts/probe_vasp_input_contract.py \
  --source-dir "${case_root}/vcneb_input/initial" \
  --trajectory "${old_workdir}/vcneb.traj" --n-images 29 \
  --workdir "${recovery_root}/initialization_probe" --symprec "${vasp_symprec}" \
  --ncores 40 --vasp-bin "${vasp_bin}" --run

for endpoint in initial final; do
  "${python}" examples/run_vcneb_vasp.py \
    --initial "${case_root}/vcneb_input/initial" --final "${case_root}/vcneb_input/final" \
    --workdir "${recovery_root}/static_${endpoint}" --n-images 29 --fmax "${fmax}" --steps 0 \
    --k "${spring}" --pressure-gpa "${pressure_gpa}" --image-workers 0 --no-climb \
    --mic --cell-interpolation linear --mapping auto --align-translation \
    --minimum-distance 1.4 --maximum-deformation 0.50 --static-only --static-endpoint "${endpoint}" \
    --vasp-isym -1 --vasp-symprec "${vasp_symprec}" --ncores 40 --vasp-bin "${vasp_bin}"
done

"${python}" examples/run_vcneb_vasp.py \
  --initial "${case_root}/vcneb_input/initial" --final "${case_root}/vcneb_input/final" \
  --workdir "${recovery_root}/vcneb_b4_to_b1_tetragonal_n29_cu17_serial" \
  --n-images 29 --fmax "${fmax}" --steps "${steps}" --k "${spring}" \
  --pressure-gpa "${pressure_gpa}" --optimizer FIRE --image-workers 0 \
  --initial-static-summary "${recovery_root}/static_initial/vasp_static_summary.json" \
  --final-static-summary "${recovery_root}/static_final/vasp_static_summary.json" \
  --resume --resume-trajectory "${old_workdir}/vcneb.traj" \
  --no-climb --mic --cell-interpolation linear --mapping auto --align-translation \
  --minimum-distance 1.4 --maximum-deformation 0.50 --vasp-isym -1 --vasp-symprec "${vasp_symprec}" \
  --ncores 40 --vasp-bin "${vasp_bin}"

"${python}" scripts/audit_vcneb_result.py \
  "${recovery_root}/vcneb_b4_to_b1_tetragonal_n29_cu17_serial" \
  --max-min-distance 1.4 --max-deformation 0.50 --fmax-target "${fmax}"
