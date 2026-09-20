#!/bin/bash
# Contract-level smoke for optional ASE backends.  This is not a material run.
# Submit only after checking sinfo/squeue and set BACKEND=lammps, cp2k, qe, or abinit.
#SBATCH -J varneb-backend
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -o backend-smoke-%j.out
#SBATCH -e backend-smoke-%j.err

set -euo pipefail

repo=${REPO:-$PWD}
python=${PYTHON:-/public/home/iai806/.conda/envs/icu/bin/python}
backend=${BACKEND:-lammps}
workdir=${WORKDIR:-validation/backend_smoke/runs/${backend}_${SLURM_JOB_ID}}

module purge
case "${backend}" in
  lammps)
    module load apps/lammps/23Jun2022-intel-21
    command_name=${COMMAND:-lmp_mpi}
    ;;
  cp2k)
    module load apps/cp2k/intelmpi/2024.1-2021-plumed
    command_name=${COMMAND:-cp2k_shell.psmp}
    ;;
  qe)
    module load apps/quantum-espresso/intelmpi/7.0
    command_name=${COMMAND:-"srun --exclusive --ntasks=1 pw.x"}
    if [[ -z "${PSEUDO_DIR:-}" ]]; then
      echo "QE requires PSEUDO_DIR containing H.SG15.PBE.UPF" >&2
      exit 2
    fi
    ;;
  abinit)
    module load compiler/intel/2021.3.0
    module load apps/abinit/intelmpi/8.6.1
    command_name=${COMMAND:-"srun --exclusive --ntasks=1 abinit"}
    if [[ -z "${PSEUDO_DIR:-}" ]]; then
      echo "ABINIT requires PSEUDO_DIR containing the reviewed pseudopotential files" >&2
      exit 2
    fi
    ;;
  *)
    echo "BACKEND must be lammps, cp2k, qe, or abinit" >&2
    exit 2
    ;;
esac

cd "${repo}"
export PYTHONPATH="${repo}:${PYTHONPATH:-}"
extra_args=()
if [[ "${backend}" == qe ]]; then
  extra_args+=(--pseudo-dir "${PSEUDO_DIR}")
fi
if [[ "${backend}" == abinit ]]; then
  extra_args+=(--pseudo-dir "${PSEUDO_DIR}" --abinit-pps "${ABINIT_PPS:-psp8}")
fi
"${python}" examples/run_optional_backend_smoke.py "${backend}" \
  --command "${command_name}" --workdir "${workdir}" "${extra_args[@]}"
