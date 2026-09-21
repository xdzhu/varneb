#!/bin/bash
# Calculator-independent ASE endpoint preparation on HF.
# The result is a CONTCAR plus endpoint_relax_summary.json; it is not a
# VC-NEB production run.  Feed both relaxed endpoints through --static-only
# and the endpoint gate before launching a path.
#SBATCH -J varneb-endpoint
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260921-varneb-material/slurm-%x-%j.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260921-varneb-material/slurm-%x-%j.err

set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 ]] || { echo "set RUN_DFT=1" >&2; exit 2; }

repo=${REPO:-/public/home/iai806/abacus/agent-runs/20260921-varneb-material/varneb_guarded_source}
python=${PYTHON:-/public/home/iai806/.conda/envs/icu/bin/python}
backend=${BACKEND:?set BACKEND=cp2k,qe,abinit,or lammps}
structure=${STRUCTURE:?set STRUCTURE endpoint path}
workdir=${WORKDIR:?set WORKDIR}
parameters=${PARAMETERS_JSON:?set PARAMETERS_JSON JSON file}
factory_kwargs=${FACTORY_KWARGS_JSON:-}
command=${COMMAND:-}
pressure=${PRESSURE_GPA:-0.0}
fmax=${FMAX:-0.10}
steps=${STEPS:-200}
maxstep=${MAXSTEP:-0.05}

module purge
case "${backend}" in
  cp2k)
    module load apps/cp2k/intelmpi/2024.1-2021-plumed
    factory=${FACTORY:-vcneb.backends:make_ase_cp2k_factory}
    command=${command:-"srun --exclusive --nodes=1 --ntasks=1 --ntasks-per-node=1 cp2k_shell.psmp"}
    ;;
  qe)
    module load apps/quantum-espresso/intelmpi/7.0
    factory=${FACTORY:-vcneb.qe:make_ase_espresso_factory}
    command=${command:-"srun --exclusive --nodes=1 --ntasks=32 --ntasks-per-node=32 pw.x"}
    ;;
  abinit)
    # Match ABINIT 8.6.1's Intel 2017 Fortran/MPI ABI.  The newer compiler
    # module and the site srun wrapper are incompatible with this binary.
    module load compiler/intel/2017.5.239
    module load mpi/intelmpi/2017.4.239
    module load apps/abinit/intelmpi/8.6.1
    factory=${FACTORY:-vcneb.backends:make_ase_abinit_factory}
    # Intel MPI Hydra bootstraps through Slurm without passing its pmi_args
    # option to ABINIT.  Keep the site PMI overrides out of this launch.
    unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS I_MPI_HYDRA_BOOTSTRAP
    command=${command:-"mpiexec.hydra -bootstrap slurm -n 32 abinit"}
    ;;
  lammps)
    module load apps/lammps/23Jun2022-intel-21
    factory=${FACTORY:-vcneb.backends:make_ase_lammps_factory}
    command=${command:-"srun --exclusive --nodes=1 --ntasks=32 --ntasks-per-node=32 lmp_mpi"}
    ;;
  *)
    echo "BACKEND must be cp2k, qe, abinit, or lammps" >&2
    exit 2
    ;;
esac

cd "${repo}"
unset PYTHONHOME PYTHONPATH
export OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
args=(scripts/relax_ase_endpoint.py
  --structure "${structure}" --workdir "${workdir}"
  --factory "${factory}" --parameters "${parameters}"
  --command "${command}" --pressure-gpa "${pressure}"
  --fmax "${fmax}" --steps "${steps}" --maxstep "${maxstep}")
if [[ -n "${factory_kwargs}" ]]; then
  args+=(--factory-kwargs "${factory_kwargs}")
fi
if [[ "${FIXED_CELL:-0}" == 1 ]]; then
  args+=(--fixed-cell)
fi
if [[ "${RESUME:-0}" == 1 ]]; then
  args+=(--resume)
fi
exec "${python}" "${args[@]}"
