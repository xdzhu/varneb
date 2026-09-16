#!/bin/bash
# Seven-image BST50 VCA-VCNEB with cached relaxed endpoints on cu17.

set -euo pipefail

repo=${REPO:-/home/zhuxd/abacus/agent-runs/20260916-varneb-v}
python=${PYTHON:-/home/zhuxd/Software/anaconda3/envs/icu/bin/python}
vasp_bin=${VASP_BIN:-/home/zhuxd/Software/src/vasp/6.3.2/bin/vasp_std}
relax_root=${RELAX_ROOT:-${repo}/batio3_vasp_vca_bst50_relax}
initial_dir=${VASP_INITIAL_DIR:-${relax_root}/vcneb_input/initial}
final_dir=${VASP_FINAL_DIR:-${relax_root}/vcneb_input/final}
initial_static=${VASP_INITIAL_STATIC:-${relax_root}/static_initial/vasp_static_summary.json}
final_static=${VASP_FINAL_STATIC:-${relax_root}/static_final/vasp_static_summary.json}
workdir=${WORKDIR:-${relax_root}/vcneb_tetragonal_to_cubic_n7_cu17_serial_db3b80f}
steps=${STEPS:-300}
fmax=${FMAX:-0.10}
spring=${SPRING:-0.20}

[[ "${RUN_DFT:-0}" == 1 ]] || { echo "Set RUN_DFT=1 to run VASP" >&2; exit 2; }
[[ "$(hostname -s)" == cu17 ]] || { echo "This runner is restricted to cu17" >&2; exit 2; }
[[ "$(nproc)" == 40 ]] || { echo "cu17 VCNEB runner requires exactly 40 visible CPU cores" >&2; exit 2; }
[[ -x "${vasp_bin}" ]] || { echo "VASP_BIN is not executable: ${vasp_bin}" >&2; exit 2; }
: "${VCNEB_GIT_REVISION:?Set the exact source commit in VCNEB_GIT_REVISION}"
for path in "${initial_dir}/CONTCAR" "${initial_dir}/INCAR" "${initial_dir}/KPOINTS" "${initial_dir}/POTCAR" "${final_dir}/CONTCAR" "${initial_static}" "${final_static}"; do
  [[ -f "${path}" ]] || { echo "missing required VCA-VCNEB input: ${path}" >&2; exit 2; }
done

export OMP_NUM_THREADS=1
export VASP_COMMAND="mpirun -np 40 ${vasp_bin}"
cd "${repo}"
command=(
  "${python}" examples/run_vcneb_vasp.py
  --initial "${initial_dir}" --final "${final_dir}" --workdir "${workdir}"
  --n-images 7 --fmax "${fmax}" --steps "${steps}" --k "${spring}" --optimizer FIRE
  --image-workers 0 --initial-static-summary "${initial_static}" --final-static-summary "${final_static}"
  --vca-virtual-symbol Ba --vca-components Ba Sr
  --no-climb --mic --cell-interpolation log_strain --mapping auto --align-translation
  --minimum-distance 1.6 --maximum-deformation 0.10
  --ncores 40 --vasp-bin "${vasp_bin}"
)
if [[ "${RESUME:-0}" == 1 ]]; then
  command+=(--resume)
fi
"${command[@]}"

"${python}" scripts/audit_vcneb_result.py "${workdir}" \
  --max-min-distance 1.6 --max-deformation 0.10 --fmax-target "${fmax}"
