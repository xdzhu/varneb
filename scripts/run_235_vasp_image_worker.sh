#!/bin/bash
# Launch exactly one VASP interior image on its assigned 28-core PBS node.
set -euo pipefail

vasp_bin=${1:?missing VASP executable}
worker_root=${2:?missing worker-root path}
image=$(basename "$PWD")
hostfile="${worker_root}/pbs_worker_hosts/${image}.host"

case "${image}" in
  01|02|03|04|05) ;;
  *) echo "This PBS worker wrapper only accepts interior image directories 01..05, got ${image}" >&2; exit 2 ;;
esac
[[ -x "${vasp_bin}" ]] || { echo "VASP executable is not available: ${vasp_bin}" >&2; exit 2; }
[[ -f "${hostfile}" ]] || { echo "missing assigned PBS host file: ${hostfile}" >&2; exit 2; }
host=$(<"${hostfile}")
[[ -n "${host}" ]] || { echo "empty assigned PBS host file: ${hostfile}" >&2; exit 2; }

exec mpirun -np 28 -hosts "${host}" "${vasp_bin}"
