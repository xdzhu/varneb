#!/bin/bash
# One bounded G2 chain, seven internal images and two cached clamped wells.
#SBATCH -J hfo2-G2-pilot
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=04:00:00
set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_JOB_NUM_NODES:-0}" == 1 \
   && "${SLURM_NTASKS:-0}" == 1 && "${SLURM_CPUS_PER_TASK:-0}" == 32 ]] || exit 2
source_root=${SOURCE_ROOT:?set fresh immutable tested archive}
seed_root=${SEED_ROOT:?set audited prepared G2 chain}
workdir=${WORKDIR:?set fresh optimization workdir}
[[ ! -e "${workdir}" ]] || { echo "refusing existing run" >&2; exit 2; }
module purge
source /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS
export I_MPI_HYDRA_BOOTSTRAP=fork
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONPATH="${source_root}${PYTHONPATH:+:${PYTHONPATH}}"
cd "${source_root}"
cell_scale=$(/public/home/iai806/.conda/envs/icu/bin/python -c \
 'import json,sys; m=json.load(open(sys.argv[1])); assert m["status"]=="G2_geometry_prepared_DFT_pending" and m["n_total_images"]==9 and m["n_internal_images"]==7 and m["n_fixed_endpoints"]==2 and m["fmax_eV_A"]==.10 and not m["climb"] and not m["physical_inputs_changed"] and not m["holdout_generated"] and m["strain"] in (0,.01); print(m["cell_scale_A"])' "${seed_root}/manifest.json")
/public/home/iai806/.conda/envs/icu/bin/python -m vcneb.material_runner \
 --initial "${seed_root}/initial.vasp" --final "${seed_root}/final.vasp" \
 --initial-chain "${seed_root}/seed.traj" --n-images 9 --workdir "${workdir}" \
 --factory examples.hfo2_fixed_input_factory:make_clamped_seed_cached_factory \
 --parameters "${seed_root}/factory_parameters.json" \
 --command "mpirun -np 32 /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus" \
 --clamped-plane-reference "${seed_root}/substrate.vasp" --clamped-allow-tilt true --cell-scale "${cell_scale}" \
 --pressure-gpa 0 --cell-mode full --optimizer FIRE --fmax .10 --steps 10 --maxstep .02 --k .2 \
 --no-climb --mapping identity --mic --no-align-cells --cell-interpolation linear \
 --minimum-distance 1.6 --maximum-deformation .25 --maximum-cell-step .02 --candidate-step-retries 4 \
 --require-continuous-periodic-lift --summary-format brief
