#!/bin/bash
#SBATCH -J hfo2-prl-end
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
set -eo pipefail
source /public/home/iai806/Software/VASP/env.sh 6.3.2
set -u

[[ "${RUN_DFT:-0}" == 1 && "${SLURM_NTASKS}" == 32 ]] || exit 2
case "${CASE_NAME:-}" in t_to_po|po_to_m) ;; *) exit 2 ;; esac
case "${ENDPOINT:-}" in initial|final) ;; *) exit 2 ;; esac

repo=/public/home/iai806/abacus/agent-runs/20260919-prl2023-hfo2-fig2a
python=/public/home/iai806/.conda/envs/icu/bin/python
case_root="${repo}/cases/${CASE_NAME}"
source_dir="${case_root}/input/${ENDPOINT}"
workdir="${case_root}/endpoint_relax/${ENDPOINT}"
export OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
export VASP_COMMAND="srun --exclusive --exact --nodes=1 --ntasks=32 --ntasks-per-node=32 --cpus-per-task=1 --mpi=pmi2 ${VASP_ROOT}/bin/vasp_std"
cd "${repo}"

for required in POSCAR INCAR KPOINTS POTCAR; do
  [[ -s "${source_dir}/${required}" ]] || { echo "missing ${source_dir}/${required}" >&2; exit 3; }
done

run_stage() {
  local target=$1
  local attempt rc resume_flag
  for attempt in 1 2 3; do
    resume_flag=""
    [[ -s "${workdir}/CONTCAR.current" ]] && resume_flag="--resume"
    set +e
    "${python}" scripts/relax_vasp_endpoint.py \
      --structure "${source_dir}/POSCAR" --source-dir "${source_dir}" --workdir "${workdir}" \
      --pressure-gpa 0 --fmax "${target}" --steps 160 --maxstep 0.03 \
      --symprec 1e-4 --ncores 32 --vasp-bin "${VASP_ROOT}/bin/vasp_std" ${resume_flag:+--resume}
    rc=$?
    set -e
    [[ ${rc} == 0 || ${rc} == 1 ]] || exit "${rc}"
    if "${python}" -c 'import json,sys; sys.exit(0 if json.load(open(sys.argv[1]))["converged"] else 1)' "${workdir}/endpoint_relax_summary.json"; then
      return 0
    fi
  done
  echo "endpoint failed to reach ${target} eV/A after three chunks" >&2
  return 1
}

run_stage 0.02
run_stage 0.005
run_stage 0.001
cp "${source_dir}/INCAR" "${source_dir}/KPOINTS" "${source_dir}/POTCAR" "${workdir}/"
sha256sum "${source_dir}/INCAR" "${source_dir}/KPOINTS" "${source_dir}/POTCAR" > "${workdir}/input.sha256"
