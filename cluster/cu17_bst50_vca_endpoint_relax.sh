#!/bin/bash
# Sequential symmetry-preserving BST50 T/C endpoint relaxation on cu17.

set -euo pipefail

repo=${REPO:-/home/zhuxd/abacus/agent-runs/20260916-varneb-v}
python=${PYTHON:-/home/zhuxd/Software/anaconda3/envs/icu/bin/python}
vasp_bin=${VASP_BIN:-/home/zhuxd/Software/src/vasp/6.3.2/bin/vasp_std}
case_root=${CASE_ROOT:-${repo}/batio3_vasp_vca_bst50_static}
relax_root=${RELAX_ROOT:-${repo}/batio3_vasp_vca_bst50_relax}

[[ "${RUN_DFT:-0}" == 1 ]] || { echo "Set RUN_DFT=1 to run VASP" >&2; exit 2; }
[[ "$(hostname -s)" == cu17 ]] || { echo "This runner is restricted to cu17" >&2; exit 2; }
[[ "$(nproc)" == 40 ]] || { echo "cu17 endpoint runner requires exactly 40 visible CPU cores" >&2; exit 2; }
[[ -x "${vasp_bin}" ]] || { echo "VASP_BIN is not executable: ${vasp_bin}" >&2; exit 2; }
: "${VCNEB_GIT_REVISION:?Set the exact source commit in VCNEB_GIT_REVISION}"

export OMP_NUM_THREADS=1
export VASP_COMMAND="mpirun -np 40 ${vasp_bin}"
cd "${repo}"

for endpoint in tetragonal cubic; do
  if [[ "${endpoint}" == tetragonal ]]; then
    structure=${case_root}/initial/CONTCAR
  else
    structure=${case_root}/final/CONTCAR
  fi
  workdir=${relax_root}/${endpoint}
  echo "START ${endpoint} $(date --iso-8601=seconds)"
  "${python}" scripts/relax_vasp_vca_endpoint.py \
    --structure "${structure}" --source-dir "${case_root}/initial" --workdir "${workdir}" \
    --virtual-symbol Ba --components Ba Sr --fmax 0.03 --steps 100 --maxstep 0.08 \
    --ncores 40 --vasp-bin "${vasp_bin}"
  echo "DONE ${endpoint} $(date --iso-8601=seconds)"
done
