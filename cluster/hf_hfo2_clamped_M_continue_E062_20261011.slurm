#!/bin/bash
# Same independent PO-to-M G2 chain; audited step-20 frame, no input tuning.
# Runtime remains the immutable E054 tested source; FIRE history is not restored.
#SBATCH -J hfo2-G2-M-E062
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=06:00:00
set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_JOB_NUM_NODES:-0}" == 1 \
   && "${SLURM_NTASKS:-0}" == 1 && "${SLURM_CPUS_PER_TASK:-0}" == 32 ]] || exit 2
source_root=/public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-G2-E054-r1/source-fixed
seed_root=${SEED_ROOT:?set audited step-20 continuation}
workdir=${WORKDIR:?set fresh workdir}
[[ ! -e "${workdir}" ]] || exit 2
module purge
source /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS
export I_MPI_HYDRA_BOOTSTRAP=fork
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONPATH="${source_root}${PYTHONPATH:+:${PYTHONPATH}}"
cd "${source_root}"
cell_scale=$(/public/home/iai806/.conda/envs/icu/bin/python -c \
 'import json,sys; m=json.load(open(sys.argv[1])); assert m["status"]=="audited_G2_geometry_continuation_prepared" and m["source_job_id"]=="28661020" and m["source_step"]==20 and m["current_frame_exact_caches"]==9 and m["new_independent_chains"]==0 and not m["physical_inputs_changed"] and not m["climb"] and not m["holdout_generated_or_read"]; print(m["cell_scale_A"])' "${seed_root}/manifest.json")
/public/home/iai806/.conda/envs/icu/bin/python -m vcneb.material_runner \
 --initial "${seed_root}/initial.vasp" --final "${seed_root}/final.vasp" \
 --initial-chain "${seed_root}/seed.traj" --n-images 9 --workdir "${workdir}" \
 --factory examples.hfo2_fixed_input_factory:make_clamped_resume_cached_factory \
 --parameters "${seed_root}/factory_parameters.json" \
 --command "mpirun -np 32 /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus" \
 --clamped-plane-reference "${seed_root}/substrate.vasp" --clamped-allow-tilt true --cell-scale "${cell_scale}" \
 --pressure-gpa 0 --cell-mode full --optimizer FIRE --fmax .10 --steps 20 --maxstep .02 --k .2 \
 --no-climb --mapping identity --mic --no-align-cells --cell-interpolation linear \
 --minimum-distance 1.6 --maximum-deformation .25 --maximum-cell-step .02 --candidate-step-retries 4 \
 --require-continuous-periodic-lift --summary-format brief
