#!/bin/bash
# Run PTO then PZT50-VCA T-to-C workflows serially on the owned cu17 node.

set -euo pipefail

repo=${REPO:-/home/zhuxd/abacus/agent-runs/20260916-varneb-v}
python=${PYTHON:-/home/zhuxd/Software/anaconda3/envs/icu/bin/python}
vasp_bin=${VASP_BIN:-/home/zhuxd/Software/src/vasp/6.3.2/bin/vasp_std}
pto_case=${PTO_CASE:-${repo}/ptio3_vasp_pbe_paw_static}
pto_relax=${PTO_RELAX:-${repo}/ptio3_vasp_relax}
pzt_case=${PZT_CASE:-${repo}/pzt50_vasp_vca_static}
pzt_relax=${PZT_RELAX:-${repo}/pzt50_vasp_vca_relax}
steps=${STEPS:-300}
fmax=${FMAX:-0.10}
spring=${SPRING:-0.20}

[[ "${RUN_DFT:-0}" == 1 ]] || { echo "Set RUN_DFT=1 to run VASP" >&2; exit 2; }
[[ "$(hostname -s)" == cu17 ]] || { echo "This runner is restricted to cu17" >&2; exit 2; }
[[ "$(nproc)" == 40 ]] || { echo "cu17 serial workflow requires exactly 40 visible CPU cores" >&2; exit 2; }
[[ -x "${vasp_bin}" ]] || { echo "VASP_BIN is not executable: ${vasp_bin}" >&2; exit 2; }
: "${VCNEB_GIT_REVISION:?Set the exact source commit in VCNEB_GIT_REVISION}"

export OMP_NUM_THREADS=1
export VASP_COMMAND="mpirun -np 40 ${vasp_bin}"
cd "${repo}"

completed() {
  [[ -f "$1" ]] && grep -q '"status": "completed"' "$1"
}

wait_or_resume_pto_relax() {
  while ! completed "${pto_relax}/tetragonal/endpoint_relax_summary.json" || ! completed "${pto_relax}/cubic/endpoint_relax_summary.json"; do
    if pgrep -f 'relax_vasp_vca_endpoint.py.*ptio3_vasp_relax' >/dev/null; then
      sleep 60
      continue
    fi
    for endpoint in tetragonal cubic; do
      summary=${pto_relax}/${endpoint}/endpoint_relax_summary.json
      completed "${summary}" && continue
      [[ "${endpoint}" == tetragonal ]] && structure=${pto_case}/initial/CONTCAR || structure=${pto_case}/final/CONTCAR
      resume=()
      [[ -d "${pto_relax}/${endpoint}" ]] && resume=(--resume)
      "${python}" scripts/relax_vasp_vca_endpoint.py \
        --structure "${structure}" --source-dir "${pto_case}/initial" --workdir "${pto_relax}/${endpoint}" \
        --fmax 0.03 --steps 100 --maxstep 0.08 --ncores 40 --vasp-bin "${vasp_bin}" "${resume[@]}"
    done
  done
}

prepare_vcneb_input() {
  local relax_root=$1 input_root=$2
  if [[ ! -d "${input_root}" ]]; then
    mkdir -p "${input_root}/initial" "${input_root}/final"
    cp "${relax_root}/tetragonal/CONTCAR" "${input_root}/initial/CONTCAR"
    cp "${relax_root}/cubic/CONTCAR" "${input_root}/final/CONTCAR"
    cp "${relax_root}/../$(basename "${relax_root}" | sed 's/_relax.*//')_vasp_pbe_paw_static/initial/INCAR" "${input_root}/initial/INCAR" 2>/dev/null || true
  fi
}

copy_vcneb_inputs() {
  local case_root=$1 relax_root=$2 input_root=${relax_root}/vcneb_input
  if [[ ! -d "${input_root}" ]]; then
    mkdir -p "${input_root}/initial" "${input_root}/final"
    cp "${relax_root}/tetragonal/CONTCAR" "${input_root}/initial/CONTCAR"
    cp "${relax_root}/cubic/CONTCAR" "${input_root}/final/CONTCAR"
    for name in INCAR KPOINTS POTCAR; do cp "${case_root}/initial/${name}" "${input_root}/initial/${name}"; done
  fi
}

run_static() {
  local case_root=$1 relax_root=$2 endpoint=$3 virtual=${4:-} components=${5:-}
  local input_root=${relax_root}/vcneb_input static_root=${relax_root}/static_${endpoint}
  completed "${static_root}/vasp_static_summary.json" && return
  local args=(--initial "${input_root}/initial" --final "${input_root}/final" --workdir "${static_root}"
    --n-images 7 --steps 0 --static-only --static-endpoint "${endpoint}" --mic --cell-interpolation log_strain
    --mapping auto --align-translation --minimum-distance 1.6 --maximum-deformation 0.10 --ncores 40 --vasp-bin "${vasp_bin}")
  [[ -n "${virtual}" ]] && args+=(--vca-virtual-symbol "${virtual}" --vca-components ${components})
  "${python}" examples/run_vcneb_vasp.py "${args[@]}"
}

run_vcneb() {
  local relax_root=$1 label=$2 virtual=${3:-} components=${4:-}
  local input_root=${relax_root}/vcneb_input workdir=${relax_root}/vcneb_tetragonal_to_cubic_n7_cu17_serial
  completed "${workdir}/vcneb_summary.json" && return
  local args=(--initial "${input_root}/initial" --final "${input_root}/final" --workdir "${workdir}"
    --n-images 7 --fmax "${fmax}" --steps "${steps}" --k "${spring}" --optimizer FIRE --image-workers 0
    --initial-static-summary "${relax_root}/static_initial/vasp_static_summary.json"
    --final-static-summary "${relax_root}/static_final/vasp_static_summary.json"
    --no-climb --mic --cell-interpolation log_strain --mapping auto --align-translation
    --minimum-distance 1.6 --maximum-deformation 0.10 --ncores 40 --vasp-bin "${vasp_bin}")
  [[ -f "${workdir}/vcneb.traj" ]] && args+=(--resume)
  [[ -n "${virtual}" ]] && args+=(--vca-virtual-symbol "${virtual}" --vca-components ${components})
  "${python}" examples/run_vcneb_vasp.py "${args[@]}"
  "${python}" scripts/audit_vcneb_result.py "${workdir}" --max-min-distance 1.6 --max-deformation 0.10 --fmax-target "${fmax}"
  echo "[DONE] ${label} VCNEB"
}

wait_or_resume_pto_relax
copy_vcneb_inputs "${pto_case}" "${pto_relax}"
run_static "${pto_case}" "${pto_relax}" initial
run_static "${pto_case}" "${pto_relax}" final
run_vcneb "${pto_relax}" PTO

if [[ ! -f "${pzt_case}/setup_manifest.json" ]]; then
  "${python}" scripts/setup_perovskite_vasp_case.py \
    --template-initial "${pto_relax}/tetragonal/CONTCAR" --template-final "${pto_relax}/cubic/CONTCAR" \
    --destination "${pzt_case}" --incar "${pto_case}/initial/INCAR" --kpoints "${pto_case}/initial/KPOINTS" \
    --template-a Pb --template-b Ti --a-symbol Pb --b-symbol Ti \
    --a-potcar /home/zhuxd/vasp/PSEUDO/PBE/Pb_d/POTCAR --b-potcar /home/zhuxd/vasp/PSEUDO/PBE/Ti_sv/POTCAR \
    --vca-b-potcar /home/zhuxd/vasp/PSEUDO/PBE/Zr_sv/POTCAR --vca-b-symbol Zr --o-potcar /home/zhuxd/vasp/PSEUDO/PBE/O/POTCAR
fi

for endpoint in tetragonal cubic; do
  completed "${pzt_relax}/${endpoint}/endpoint_relax_summary.json" && continue
  [[ "${endpoint}" == tetragonal ]] && structure=${pzt_case}/initial/CONTCAR || structure=${pzt_case}/final/CONTCAR
  resume=()
  [[ -d "${pzt_relax}/${endpoint}" ]] && resume=(--resume)
  "${python}" scripts/relax_vasp_vca_endpoint.py \
    --structure "${structure}" --source-dir "${pzt_case}/initial" --workdir "${pzt_relax}/${endpoint}" \
    --virtual-symbol Ti --components Ti Zr --fmax 0.03 --steps 100 --maxstep 0.08 \
    --ncores 40 --vasp-bin "${vasp_bin}" "${resume[@]}"
done
copy_vcneb_inputs "${pzt_case}" "${pzt_relax}"
run_static "${pzt_case}" "${pzt_relax}" initial Ti "Ti Zr"
run_static "${pzt_case}" "${pzt_relax}" final Ti "Ti Zr"
run_vcneb "${pzt_relax}" PZT50 Ti "Ti Zr"
