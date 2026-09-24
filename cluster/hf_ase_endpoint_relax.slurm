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
stress_kbar=${STRESS_KBAR:-1.0}
abinit_mpi_ranks=${ABINIT_MPI_RANKS:-8}
qe_mpi_ranks=${QE_MPI_RANKS:-1}
cp2k_mpi_ranks=${CP2K_MPI_RANKS:-16}

module purge
case "${backend}" in
  cp2k)
    module load apps/cp2k/intelmpi/2024.1-2021-plumed
    factory=${FACTORY:-vcneb.backends:make_ase_cp2k_factory}
    # cp2k_shell is one logical ASE service, but its PSMP worker may and should
    # be launched as one MPI world.  On HF, direct Intel-MPI launch is both
    # correct and measured faster than the single-rank shell (16 ranks: 2m23s;
    # 8 ranks: 3m44s for the same GaN 800-Ry static calculation).
    command=${command:-"mpirun -np ${cp2k_mpi_ranks} cp2k_shell.psmp"}
    ;;
  qe)
    module load apps/quantum-espresso/intelmpi/7.0
    factory=${FACTORY:-vcneb.qe:make_ase_espresso_factory}
    # HF's QE 7.0 binary reaches MPI_Init reliably through the site Slurm
    # launcher only as a single-rank step.  A 32-rank srun step fails before
    # reading the input (PMI/proxy socket startup), so keep the proven launch
    # contract explicit instead of mistaking that infrastructure failure for
    # a pseudopotential or SCF problem.
    command=${command:-"srun --exclusive --nodes=1 --ntasks=${qe_mpi_ranks} --ntasks-per-node=${qe_mpi_ranks} pw.x"}
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
    command=${command:-"mpiexec.hydra -bootstrap slurm -n ${abinit_mpi_ranks} abinit"}
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
  --fmax "${fmax}" --stress-kbar "${stress_kbar}"
  --steps "${steps}" --maxstep "${maxstep}")
if [[ -n "${factory_kwargs}" ]]; then
  args+=(--factory-kwargs "${factory_kwargs}")
fi
if [[ "${FIXED_CELL:-0}" == 1 ]]; then
  # A fixed-cell smoke test cannot satisfy a variable-cell stress target.
  # Rebuild the argument vector without the stress option and its value.
  args=(scripts/relax_ase_endpoint.py
    --structure "${structure}" --workdir "${workdir}"
    --factory "${factory}" --parameters "${parameters}"
    --command "${command}" --pressure-gpa "${pressure}"
    --fmax "${fmax}" --steps "${steps}" --maxstep "${maxstep}")
  if [[ -n "${factory_kwargs}" ]]; then
    args+=(--factory-kwargs "${factory_kwargs}")
  fi
  args+=(--fixed-cell)
fi
if [[ "${RESUME:-0}" == 1 ]]; then
  args+=(--resume)
fi
exec "${python}" "${args[@]}"
