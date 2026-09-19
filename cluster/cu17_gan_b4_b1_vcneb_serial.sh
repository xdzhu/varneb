#!/bin/bash
# Serial VASP workflow for the Qian GaN B4(wurtzite) -> B1(rocksalt) VCNEB.
#
# The endpoint BFGS relaxations and every interior image use all 40 cu17
# cores.  Interior images are deliberately serial (image-workers=0); the two
# endpoint static SCFs are completed once and then cached by the manager.

set -euo pipefail

repo=${REPO:-/home/zhuxd/abacus/agent-runs/20260916-varneb-v}
python=${PYTHON:-/home/zhuxd/Software/anaconda3/envs/icu/bin/python}
vasp_bin=${VASP_BIN:-/home/zhuxd/Software/src/vasp/6.3.2/bin/vasp_std}
smoke_template=${GAN_SMOKE_TEMPLATE:-/home/zhuxd/abacus/agent-runs/20260912-vcneb-p0/validation/vasp_gan_smoke/template}
case_root=${CASE_ROOT:-${repo}/gan_b4_b1_vasp_pbe_paw_qian}
pressure_gpa=${PRESSURE_GPA:-45.7}
steps=${STEPS:-300}
fmax=${FMAX:-0.10}
spring=${SPRING:-0.20}

[[ "${RUN_DFT:-0}" == 1 ]] || { echo "Set RUN_DFT=1 to run VASP on cu17" >&2; exit 2; }
[[ "$(hostname -s)" == cu17 ]] || { echo "This runner is restricted to cu17" >&2; exit 2; }
visible_cores=$(env -u OMP_NUM_THREADS -u OMP_THREAD_LIMIT nproc)
[[ "${visible_cores}" == 40 ]] || { echo "cu17 serial workflow requires exactly 40 visible CPU cores" >&2; exit 2; }
[[ -x "${vasp_bin}" ]] || { echo "VASP_BIN is not executable: ${vasp_bin}" >&2; exit 2; }
: "${VCNEB_GIT_REVISION:?Set the exact source commit in VCNEB_GIT_REVISION}"
for path in "${smoke_template}/POSCAR" "${smoke_template}/INCAR" "${smoke_template}/KPOINTS" "${smoke_template}/POTCAR"; do
  [[ -f "${path}" ]] || { echo "missing GaN template input: ${path}" >&2; exit 2; }
done
[[ ! -e "${case_root}" ]] || { echo "case root already exists; use a fresh CASE_ROOT or an explicit recovery workflow: ${case_root}" >&2; exit 2; }

export OMP_NUM_THREADS=1
export VASP_COMMAND="mpirun -np 40 ${vasp_bin}"
cd "${repo}"

"${python}" scripts/setup_gan_b4_b1_case.py \
  --b4-poscar "${smoke_template}/POSCAR" --source-dir "${smoke_template}" --destination "${case_root}"

for endpoint in B4 B1; do
  "${python}" scripts/relax_vasp_endpoint.py \
    --structure "${case_root}/template/${endpoint}.POSCAR" --source-dir "${case_root}/template" \
    --workdir "${case_root}/endpoint_relax_${endpoint}" --pressure-gpa "${pressure_gpa}" \
    --fmax 0.02 --steps 120 --maxstep 0.05 --ncores 40 --vasp-bin "${vasp_bin}"
done

mkdir -p "${case_root}/vcneb_input/initial" "${case_root}/vcneb_input/final"
cp "${case_root}/endpoint_relax_B4/CONTCAR" "${case_root}/vcneb_input/initial/CONTCAR"
cp "${case_root}/endpoint_relax_B1/CONTCAR" "${case_root}/vcneb_input/final/CONTCAR"
cp "${case_root}/template/INCAR" "${case_root}/template/KPOINTS" "${case_root}/template/POTCAR" "${case_root}/vcneb_input/initial/"

for endpoint in initial final; do
  "${python}" examples/run_vcneb_vasp.py \
    --initial "${case_root}/vcneb_input/initial" --final "${case_root}/vcneb_input/final" \
    --workdir "${case_root}/static_${endpoint}" --n-images 29 --fmax "${fmax}" --steps 0 --k "${spring}" \
    --pressure-gpa "${pressure_gpa}" --image-workers 0 --no-climb --mic --cell-interpolation linear \
    --mapping auto --align-translation --minimum-distance 1.4 --maximum-deformation 0.50 \
    --static-only --static-endpoint "${endpoint}" --vasp-isym -1 --ncores 40 --vasp-bin "${vasp_bin}"
done

"${python}" examples/run_vcneb_vasp.py \
  --initial "${case_root}/vcneb_input/initial" --final "${case_root}/vcneb_input/final" \
  --workdir "${case_root}/vcneb_b4_to_b1_tetragonal_n29_cu17_serial" \
  --n-images 29 --fmax "${fmax}" --steps "${steps}" --k "${spring}" --pressure-gpa "${pressure_gpa}" \
  --optimizer FIRE --image-workers 0 \
  --initial-static-summary "${case_root}/static_initial/vasp_static_summary.json" \
  --final-static-summary "${case_root}/static_final/vasp_static_summary.json" \
  --no-climb --mic --cell-interpolation linear --mapping auto --align-translation \
  --minimum-distance 1.4 --maximum-deformation 0.50 --vasp-isym -1 --ncores 40 --vasp-bin "${vasp_bin}"

"${python}" scripts/audit_vcneb_result.py "${case_root}/vcneb_b4_to_b1_tetragonal_n29_cu17_serial" \
  --max-min-distance 1.4 --max-deformation 0.50 --fmax-target "${fmax}"
