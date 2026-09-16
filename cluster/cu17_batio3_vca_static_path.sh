#!/bin/bash
# Evaluate a prepared seven-image VASP VCA path serially on the owned cu17 node.

set -euo pipefail

root=${VCA_PATH_ROOT:?Set VCA_PATH_ROOT to a prepared 00..06 path}
vasp_bin=${VASP_BIN:-/home/zhuxd/Software/src/vasp/6.3.2/bin/vasp_std}

[[ "${RUN_DFT:-0}" == 1 ]] || { echo "Set RUN_DFT=1 to run VASP" >&2; exit 2; }
[[ "$(hostname -s)" == cu17 ]] || { echo "This runner is restricted to cu17" >&2; exit 2; }
[[ "$(nproc)" == 40 ]] || { echo "cu17 path runner requires exactly 40 visible CPU cores" >&2; exit 2; }
[[ -x "${vasp_bin}" ]] || { echo "VASP_BIN is not executable: ${vasp_bin}" >&2; exit 2; }

export OMP_NUM_THREADS=1
for index in 00 01 02 03 04 05 06; do
  image_dir=${root}/${index}
  for name in INCAR KPOINTS POSCAR POTCAR; do
    [[ -f "${image_dir}/${name}" ]] || { echo "missing ${image_dir}/${name}" >&2; exit 2; }
  done
  [[ ! -e "${image_dir}/OUTCAR" ]] || { echo "refusing to overwrite ${image_dir}/OUTCAR" >&2; exit 2; }
  echo "START ${index} $(date +%T)"
  cd "${image_dir}"
  mpirun -np 40 "${vasp_bin}" > vasp.out 2>&1
  grep "1 F=" OSZICAR | tail -1
  echo "DONE ${index} $(date +%T)"
done
