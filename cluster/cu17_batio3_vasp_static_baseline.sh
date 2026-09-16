#!/bin/bash
# Direct fixed-endpoint VASP static baseline for the owned 40-core cu17 node.

set -euo pipefail

repo=${REPO:-/home/zhuxd/abacus/agent-runs/20260916-varneb-v}
python=${PYTHON:-/home/zhuxd/Software/anaconda3/envs/icu/bin/python}
vasp_bin=${VASP_BIN:-/home/zhuxd/Software/src/vasp/6.3.2/bin/vasp_std}
initial_dir=${VASP_INITIAL_DIR:-${repo}/batio3_vasp_pbe_paw_static/initial}
final_dir=${VASP_FINAL_DIR:-${repo}/batio3_vasp_pbe_paw_static/final}
endpoint_identity_gate=${ENDPOINT_IDENTITY_GATE:-${repo}/batio3_vasp_pbe_paw_static/endpoint_identity_gate.json}
static_endpoint=${STATIC_ENDPOINT:-initial}
workdir=${WORKDIR:-batio3_vasp_pbe_paw_static/static_${static_endpoint}_ecut600}

[[ "${RUN_DFT:-0}" == 1 ]] || { echo "Set RUN_DFT=1 to run VASP on cu17" >&2; exit 2; }
[[ "$(hostname -s)" == cu17 ]] || { echo "This runner is restricted to cu17" >&2; exit 2; }
[[ "$(nproc)" == 40 ]] || { echo "cu17 static runner requires exactly 40 visible CPU cores" >&2; exit 2; }
[[ "${static_endpoint}" == initial || "${static_endpoint}" == final ]] || { echo "STATIC_ENDPOINT must be initial or final" >&2; exit 2; }
[[ -x "${vasp_bin}" ]] || { echo "VASP_BIN is not executable: ${vasp_bin}" >&2; exit 2; }
: "${VCNEB_GIT_REVISION:?Set the exact source commit in VCNEB_GIT_REVISION}"
for path in "${initial_dir}/CONTCAR" "${initial_dir}/INCAR" "${initial_dir}/KPOINTS" "${initial_dir}/POTCAR" "${final_dir}/CONTCAR" "${endpoint_identity_gate}"; do
  [[ -f "${path}" ]] || { echo "missing required VASP input: ${path}" >&2; exit 2; }
done

export OMP_NUM_THREADS=1
export VASP_COMMAND="mpirun -np 40 ${vasp_bin}"
cd "${repo}"
"${python}" scripts/validate_endpoint_identity_gate.py --gate "${endpoint_identity_gate}"
exec "${python}" examples/run_vcneb_vasp.py \
  --initial "${initial_dir}" --final "${final_dir}" --workdir "${workdir}" \
  --n-images 7 --fmax 0.10 --steps 0 --k 0.2 --image-workers 0 \
  --no-climb --mic --cell-interpolation log_strain --mapping auto --align-translation \
  --minimum-distance 1.6 --maximum-deformation 0.10 --static-only \
  --static-endpoint "${static_endpoint}"
