#!/bin/bash
# Direct, single-node BTO VASP VCNEB production runner for 235/cu17.
#
# Seven total images are used: fixed cached endpoints 00/06 plus five
# interior images.  `--image-workers 0` is intentional: it evaluates images
# 01..05 serially, while every individual VASP calculation uses all 40 cores
# on cu17.  Do not run this on a shared or scheduled node.

set -euo pipefail

repo=${REPO:-/home/zhuxd/abacus/agent-runs/20260916-varneb-v}
python=${PYTHON:-/home/zhuxd/Software/anaconda3/envs/icu/bin/python}
vasp_bin=${VASP_BIN:-/home/zhuxd/Software/src/vasp/6.3.2/bin/vasp_std}
initial_dir=${VASP_INITIAL_DIR:-${repo}/batio3_vasp_pbe_paw_static/initial}
final_dir=${VASP_FINAL_DIR:-${repo}/batio3_vasp_pbe_paw_static/final}
endpoint_identity_gate=${ENDPOINT_IDENTITY_GATE:-${repo}/batio3_vasp_pbe_paw_static/endpoint_identity_gate.json}
initial_static=${VASP_INITIAL_STATIC:-${repo}/batio3_vasp_pbe_paw_static/static_initial_ecut600/vasp_static_summary.json}
final_static=${VASP_FINAL_STATIC:-${repo}/batio3_vasp_pbe_paw_static/static_final_ecut600/vasp_static_summary.json}
steps=${STEPS:-300}
fmax=${FMAX:-0.10}
spring=${SPRING:-0.2}
workdir=${WORKDIR:-batio3_vasp_pbe_paw_static/vcneb_tetragonal_to_cubic_n7_cu17_serial}

[[ "${RUN_DFT:-0}" == 1 ]] || { echo "Set RUN_DFT=1 to run VASP on cu17" >&2; exit 2; }
[[ "$(hostname -s)" == cu17 ]] || { echo "This runner is restricted to cu17" >&2; exit 2; }
[[ "$(nproc)" == 40 ]] || { echo "cu17 serial runner requires exactly 40 visible CPU cores" >&2; exit 2; }
[[ -x "${vasp_bin}" ]] || { echo "VASP_BIN is not executable: ${vasp_bin}" >&2; exit 2; }
: "${VCNEB_GIT_REVISION:?Set the exact source commit in VCNEB_GIT_REVISION}"
for path in "${initial_dir}/CONTCAR" "${initial_dir}/INCAR" "${initial_dir}/KPOINTS" "${initial_dir}/POTCAR" "${final_dir}/CONTCAR" "${endpoint_identity_gate}" "${initial_static}" "${final_static}"; do
  [[ -f "${path}" ]] || { echo "missing required VASP input: ${path}" >&2; exit 2; }
done

export OMP_NUM_THREADS=1
export VASP_COMMAND="mpirun -np 40 ${vasp_bin}"
cd "${repo}"
"${python}" scripts/validate_endpoint_identity_gate.py --gate "${endpoint_identity_gate}"
"${python}" scripts/validate_vasp_static_baseline_gate.py \
  --summary "${initial_static}" --endpoint-gate "${endpoint_identity_gate}" \
  --initial-dir "${initial_dir}" --git-revision "${VCNEB_GIT_REVISION}"

command=(
  "${python}" examples/run_vcneb_vasp.py
  --initial "${initial_dir}" --final "${final_dir}" --workdir "${workdir}"
  --n-images 7 --fmax "${fmax}" --steps "${steps}" --k "${spring}" --optimizer FIRE
  --image-workers 0 --initial-static-summary "${initial_static}" --final-static-summary "${final_static}"
  --no-climb --mic --cell-interpolation log_strain --mapping auto --align-translation
  --minimum-distance 1.6 --maximum-deformation 0.10
)
if [[ "${VALIDATE_ONLY:-0}" == 1 ]]; then
  command+=(--validate-only)
fi
"${command[@]}"

if [[ "${VALIDATE_ONLY:-0}" != 1 ]]; then
  "${python}" scripts/audit_vcneb_result.py "${workdir}" \
    --max-min-distance 1.6 --max-deformation 0.10 --fmax-target "${fmax}"
fi
