#!/bin/bash
# Generic ASE-factory material VCNEB production template.
# It uses three 97-core nodes: nine 32-MPI image workers in waves plus
# manager headroom. For BTO set IMAGE_WORKERS=5; for GaN set IMAGE_WORKERS=9.
# This template is intentionally input-driven and refuses to guess potentials.
# Production FIRE steps use a 0.05 relative cell trust radius with bounded
# geometry backtracking, preventing a single image from changing volume by 10%
# while the neighboring images remain on the old path.
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
endpoint_summary=${ENDPOINT_STATIC_SUMMARY:-}
resume_snapshot=${RESUME_SNAPSHOT:-}
n_images=${N_IMAGES:-7}
image_workers=${IMAGE_WORKERS:-5}
image_mpi=${IMAGE_MPI:-32}
qe_mpi_ranks=${QE_MPI_RANKS:-1}
cp2k_mpi_ranks=${CP2K_MPI_RANKS:-16}
steps=${STEPS:-300}
fmax=${FMAX:-0.10}
pressure_gpa=${PRESSURE_GPA:-0.0}
spring=${SPRING:-0.20}
optimizer=${OPTIMIZER:-FIRE}
command=${COMMAND:-}
maximum_cell_step=${MAXIMUM_CELL_STEP:-0.05}
candidate_step_retries=${CANDIDATE_STEP_RETRIES:-8}
align_translation=${ALIGN_TRANSLATION:-1}
mapping=${MAPPING:-auto}
align_cells=${ALIGN_CELLS:-1}
cell_interpolation=${CELL_INTERPOLATION:-log_strain}

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
    # This HF QE build has a verified single-rank Slurm launch contract.  Its
    # multi-rank PMI startup fails in MPI_Init_thread before pw.x reads input.
    # Parallelism therefore comes from independent image workers, not ranks
    # inside one image, unless QE_MPI_RANKS is deliberately revalidated.
    command=${command:-"srun --exclusive --nodes=1 --ntasks=${qe_mpi_ranks} --ntasks-per-node=${qe_mpi_ranks} pw.x"}
    ;;
  cp2k)
    module load apps/cp2k/intelmpi/2024.1-2021-plumed
    factory=${FACTORY:-vcneb.backends:make_ase_cp2k_factory}
    # One ASE image owns one persistent cp2k_shell MPI world.  Direct mpirun
    # avoids the earlier nested-srun/protocol failures while retaining CP2K's
    # internal parallelism; 16 ranks is the validated HF GaN default.
    command=${command:-"mpirun -np ${cp2k_mpi_ranks} cp2k_shell.psmp"}
    ;;
  abinit)
    # The HF ABINIT 8.6.1 executable is linked against the Intel 2017
    # Fortran/MPI ABI.  Loading a newer compiler or the site srun wrapper
    # leaves ranks in PMI/socket startup (or injects the unsupported
    # ``pmi_args`` option), so keep the runtime and launcher matched here.
    module load compiler/intel/2017.5.239
    module load mpi/intelmpi/2017.4.239
    module load apps/abinit/intelmpi/8.6.1
    factory=${FACTORY:-vcneb.backends:make_ase_abinit_factory}
    # Use Intel MPI's Hydra launcher explicitly.  Do not inherit the HF
    # site's external-launcher/PMI overrides; Hydra bootstraps through Slurm
    # without passing Slurm's pmi_args to ABINIT.
    unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS I_MPI_HYDRA_BOOTSTRAP
    command=${command:-"mpiexec.hydra -bootstrap slurm -n ${image_mpi} abinit"}
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
# The gate report is part of the production artifact.  Create the isolated
# work directory before validation so a passed gate cannot fail merely because
# its output parent is absent.
mkdir -p "${workdir}"
if [[ "${STATIC_ONLY:-0}" != 1 && "${REQUIRE_ENDPOINT_STATIC_GATE:-1}" == 1 ]]; then
  [[ -n "${endpoint_summary}" ]] || {
    echo "ENDPOINT_STATIC_SUMMARY is required before production VC-NEB" >&2
    exit 2
  }
  if [[ "${endpoint_summary}" != /* ]]; then
    endpoint_summary="${repo}/${endpoint_summary}"
  fi
  "${python}" scripts/validate_ase_static_gate.py "${endpoint_summary}" \
    --fmax "${ENDPOINT_FMAX:-0.10}" \
    --stress-kbar "${ENDPOINT_STRESS_KBAR:-1.0}" \
    --pressure-gpa "${pressure_gpa}" \
    --output "${workdir}/endpoint_static_gate.json"
fi
args=(examples/run_vcneb_ase.py
  --initial "${initial}" --final "${final}" --workdir "${workdir}"
  --n-images "${n_images}" --fmax "${fmax}" --steps "${steps}" --k "${spring}"
  --pressure-gpa "${pressure_gpa}"
  --optimizer "${optimizer}" --image-workers "${image_workers}"
  --candidate-step-retries "${candidate_step_retries}" --maximum-cell-step "${maximum_cell_step}"
  --command "${command}" --factory "${factory}" --parameters "${parameters}"
  --mic --cell-interpolation "${cell_interpolation}" --mapping "${mapping}"
  --maximum-deformation "${MAXIMUM_DEFORMATION:-0.50}")
if [[ "${align_translation}" == 1 ]]; then
  args+=(--align-translation)
elif [[ "${align_translation}" != 0 ]]; then
  echo "ALIGN_TRANSLATION must be 0 or 1" >&2
  exit 2
fi
if [[ "${mapping}" != auto && "${mapping}" != identity ]]; then
  echo "MAPPING must be auto or identity" >&2
  exit 2
fi
if [[ "${align_cells}" != 0 && "${align_cells}" != 1 ]]; then
  echo "ALIGN_CELLS must be 0 or 1" >&2
  exit 2
fi
if [[ "${align_cells}" == 0 ]]; then
  args+=(--no-align-cells)
fi
if [[ "${cell_interpolation}" != linear && "${cell_interpolation}" != log_strain ]]; then
  echo "CELL_INTERPOLATION must be linear or log_strain" >&2
  exit 2
fi
if [[ -n "${endpoint_summary}" && "${STATIC_ONLY:-0}" != 1 ]]; then
  args+=(--endpoint-static-summary "${endpoint_summary}")
fi
if [[ -n "${MINIMUM_DISTANCE:-}" ]]; then
  args+=(--minimum-distance "${MINIMUM_DISTANCE}")
fi
if [[ -n "${MINIMUM_ENDPOINT_SEPARATION:-}" ]]; then
  args+=(--minimum-endpoint-separation "${MINIMUM_ENDPOINT_SEPARATION}")
fi
if [[ -n "${factory_kwargs}" ]]; then
  args+=(--factory-kwargs "${factory_kwargs}")
fi
if [[ -n "${resume_snapshot}" ]]; then
  args+=(--resume-snapshot "${resume_snapshot}")
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
