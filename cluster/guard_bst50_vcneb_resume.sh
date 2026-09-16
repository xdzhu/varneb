#!/bin/bash
# Resume the cu17 BST50 VCA-VCNEB once if its attached launcher disappears.

set -euo pipefail

watched_pid=${1:?pass the active run_vcneb_vasp.py PID}
repo=${REPO:-/home/zhuxd/abacus/agent-runs/20260916-varneb-v}
root=${WORKDIR:-${repo}/batio3_vasp_vca_bst50_relax/vcneb_tetragonal_to_cubic_n7_cu17_serial_db3b80f}
runner=${RUNNER:-${repo}/cluster/cu17_bst50_vca_vcneb.sh}

while kill -0 "${watched_pid}" 2>/dev/null; do
  sleep 30
done
if [[ ! -f "${root}/vcneb_summary.json" ]]; then
  RUN_DFT=1 RESUME=1 VCNEB_GIT_REVISION=${VCNEB_GIT_REVISION:-c5a4660} bash "${runner}"
fi
