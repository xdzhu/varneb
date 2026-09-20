#!/bin/bash
# Generic ASE-factory material VCNEB production template.
# It uses three 97-core nodes: nine 32-MPI image workers in waves plus
# manager headroom. For BTO set IMAGE_WORKERS=5; for GaN set IMAGE_WORKERS=9.
# This template is intentionally input-driven and refuses to guess potentials.
#SBATCH -J varneb-material-ase
#SBATCH -p hfacnormal01
#SBATCH -N 3
#SBATCH --ntasks=291
#SBATCH --ntasks-per-node=97
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH -o /public/home/iai806/abacus/agent-runs/20260921-varneb-material/slurm-%x-%j.out
#SBATCH -e /public/home/iai806/abacus/agent-runs/20260921-varneb-material/slurm-%x-%j.err

set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 ]] || { echo "set RUN_DFT=1 after endpoint/preflight gates" >&2; exit 2; }

repo=${REPO:-/public/home/iai806/abacus/agent-runs/20260921-varneb-material}
python=${PYTHON:-/public/home/iai806/.conda/envs/icu/bin/python}
backend=${BACKEND:?set BACKEND=qe,cp2k,abinit,or lammps}
initial=${INITIAL:?set INITIAL endpoint path}
final=${FINAL:?set FINAL endpoint path}
workdir=${WORKDIR:?set WORKDIR}
parameters=${PARAMETERS_JSON:?set PARAMETERS_JSON JSON file}
factory_kwargs=${FACTORY_KWARGS_JSON:-}
n_images=${N_IMAGES:-7}
image_workers=${IMAGE_WORKERS:-5}
image_mpi=${IMAGE_MPI:-32}
steps=${STEPS:-300}
fmax=${FMAX:-0.10}
spring=${SPRING:-0.20}
optimizer=${OPTIMIZER:-FIRE}
command=${COMMAND:-}

if (( image_workers < 1 || image_workers * image_mpi > 288 )); then
  echo "IMAGE_WORKERS * IMAGE_MPI must be between 1 and 288" >&2
  exit 2
fi
if (( SLURM_NTASKS < image_workers * image_mpi )); then
  echo "allocation is smaller than requested image worker ranks" >&2
  exit 2
fi

module purge
case "${backend}" in
  qe)
    module load apps/quantum-espresso/intelmpi/7.0
    factory=${FACTORY:-vcneb.qe:make_ase_espresso_factory}
    command=${command:-"srun --exclusive --nodes=1 --ntasks=${image_mpi} --ntasks-per-node=${image_mpi} --mpi=pmi2 pw.x"}
    ;;
  cp2k)
    module load apps/cp2k/intelmpi/2024.1-2021-plumed
    factory=${FACTORY:-vcneb.backends:make_ase_cp2k_factory}
    command=${command:-"srun --exclusive --nodes=1 --ntasks=${image_mpi} --ntasks-per-node=${image_mpi} cp2k_shell.psmp"}
    ;;
  abinit)
    module load compiler/intel/2021.3.0
    module load apps/abinit/intelmpi/8.6.1
    factory=${FACTORY:-vcneb.backends:make_ase_abinit_factory}
    command=${command:-"srun --exclusive --nodes=1 --ntasks=${image_mpi} --ntasks-per-node=${image_mpi} abinit"}
    ;;
  lammps)
    module load apps/lammps/23Jun2022-intel-21
    factory=${FACTORY:-vcneb.backends:make_ase_lammps_factory}
    command=${command:-"srun --exclusive --nodes=1 --ntasks=${image_mpi} --ntasks-per-node=${image_mpi} lmp_mpi"}
    ;;
  *)
    echo "BACKEND must be qe, cp2k, abinit, or lammps" >&2
    exit 2
    ;;
esac

cd "${repo}"
# Some HF modulefiles prepend the system Python 3.7 stdlib to PYTHONPATH.
# The reproducible VARNEB environment is the explicit interpreter above;
# leaving that path active mixes incompatible importlib implementations.
unset PYTHONHOME PYTHONPATH
export OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
args=(examples/run_vcneb_ase.py
  --initial "${initial}" --final "${final}" --workdir "${workdir}"
  --n-images "${n_images}" --fmax "${fmax}" --steps "${steps}" --k "${spring}"
  --optimizer "${optimizer}" --image-workers "${image_workers}"
  --command "${command}" --factory "${factory}" --parameters "${parameters}"
  --mic --cell-interpolation log_strain --mapping auto --align-translation
  --maximum-deformation "${MAXIMUM_DEFORMATION:-0.50}")
if [[ -n "${MINIMUM_DISTANCE:-}" ]]; then
  args+=(--minimum-distance "${MINIMUM_DISTANCE}")
fi
if [[ -n "${factory_kwargs}" ]]; then
  args+=(--factory-kwargs "${factory_kwargs}")
fi
if [[ "${NO_CLIMB:-1}" == 1 ]]; then
  args+=(--no-climb)
fi
if [[ "${VALIDATE_ONLY:-0}" == 1 ]]; then
  args+=(--validate-only)
fi
if [[ "${STATIC_ONLY:-0}" == 1 ]]; then
  args+=(--static-only --static-endpoints "${STATIC_ENDPOINTS:-both}")
fi
exec "${python}" "${args[@]}"
