#!/bin/bash
# Bounded ordinary-chain pilot for explicitly prepared/audited HfO2 cases.
#SBATCH -J hfo2-gap-refine
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=04:00:00
set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_JOB_NUM_NODES:-0}" == 1 \
   && "${SLURM_NTASKS:-0}" == 1 && "${SLURM_CPUS_PER_TASK:-0}" == 32 ]] || exit 2
source_root=${SOURCE_ROOT:?set immutable archive}
seed_root=${SEED_ROOT:?set audited prepared chain}
workdir=${WORKDIR:?set fresh optimization workdir}
[[ ! -e "${workdir}" ]] || { echo "refusing existing run" >&2; exit 2; }
module purge
source /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS
export I_MPI_HYDRA_BOOTSTRAP=fork
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONPATH="${source_root}${PYTHONPATH:+:${PYTHONPATH}}"
cd "${source_root}"
n_images=$(/public/home/iai806/.conda/envs/icu/bin/python -c \
  'import json,sys; m=json.load(open(sys.argv[1])); n=m["n_total_images"]; assert n in (9,10) and m["fmax_eV_A"]==0.10 and not m["climb"] and not m["physical_inputs_changed"]; print(n)' "${seed_root}/manifest.json")
/public/home/iai806/.conda/envs/icu/bin/python -m vcneb.material_runner \
  --initial "${seed_root}/initial.vasp" --final "${seed_root}/final.vasp" \
  --initial-chain "${seed_root}/seed.traj" --n-images "${n_images}" --workdir "${workdir}" \
  --factory examples.hfo2_fixed_input_factory:make_seed_cached_factory \
  --parameters "${seed_root}/factory_parameters.json" \
  --command "mpirun -np 32 /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus" \
  --pressure-gpa 0 --cell-mode full --optimizer FIRE --fmax 0.10 --steps 10 \
  --maxstep 0.02 --k 0.2 --no-climb --mapping identity --mic \
  --no-align-cells --cell-interpolation linear --minimum-distance 1.6 \
  --maximum-deformation 0.25 --maximum-cell-step 0.02 --candidate-step-retries 4 \
  --summary-format brief
