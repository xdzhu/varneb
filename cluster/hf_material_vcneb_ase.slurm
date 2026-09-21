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
steps=${STEPS:-300}
fmax=${FMAX:-0.10}
spring=${SPRING:-0.20}
optimizer=${OPTIMIZER:-FIRE}
command=${COMMAND:-}
maximum_cell_step=${MAXIMUM_CELL_STEP:-0.05}
candidate_step_retries=${CANDIDATE_STEP_RETRIES:-8}

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
    command=${command:-"srun --exclusive --nodes=1 --ntasks=${image_mpi} --ntasks-per-node=${image_mpi} pw.x"}
    ;;
  cp2k)
    module load apps/cp2k/intelmpi/2024.1-2021-plumed
    factory=${FACTORY:-vcneb.backends:make_ase_cp2k_factory}
    # ASE's CP2K calculator speaks to one persistent cp2k_shell process per
    # image.  Do not launch the shell as a 32-rank MPI application: the shell
    # handshake is single-rank even when the CP2K allocation has more cores.
    command=${command:-"srun --exclusive --nodes=1 --ntasks=1 --ntasks-per-node=1 cp2k_shell.psmp"}
    ;;
  abinit)
    module load compiler/intel/2021.3.0
    module load apps/abinit/intelmpi/8.6.1
    factory=${FACTORY:-vcneb.backends:make_ase_abinit_factory}
    # HF's ABINIT 8.6.1 binary is launched by the Slurm-aware mpiexec
    # wrapper.  Calling it through srun injects PMI arguments that ABINIT
    # misinterprets and leaves all ranks waiting in the MPI socket layer.
    command=${command:-"mpiexec -n ${image_mpi} abinit"}
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
    --stress-kbar "${ENDPOINT_STRESS_KBAR:-0.10}" \
    --output "${workdir}/endpoint_static_gate.json"
fi
args=(examples/run_vcneb_ase.py
  --initial "${initial}" --final "${final}" --workdir "${workdir}"
  --n-images "${n_images}" --fmax "${fmax}" --steps "${steps}" --k "${spring}"
  --optimizer "${optimizer}" --image-workers "${image_workers}"
  --candidate-step-retries "${candidate_step_retries}" --maximum-cell-step "${maximum_cell_step}"
  --command "${command}" --factory "${factory}" --parameters "${parameters}"
  --mic --cell-interpolation log_strain --mapping auto --align-translation
  --maximum-deformation "${MAXIMUM_DEFORMATION:-0.50}")
if [[ -n "${MINIMUM_DISTANCE:-}" ]]; then
  args+=(--minimum-distance "${MINIMUM_DISTANCE}")
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
