#!/bin/bash
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=18:00:00
set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 && "${CHANNEL_ID:-x}" =~ ^[01]$ ]] || exit 2
root=/public/home/iai806/abacus/agent-runs/20261019-varneb-preserving-strain-E275-r1
runtime=/public/home/iai806/abacus/agent-runs/20261018-varneb-physical-E268-r1
if [[ "${CHANNEL_ID}" == 0 ]]; then label=strain_0050; else label=strain_0100; fi
seed="${root}/bundle/seeds/${label}"
work="${root}/paths/${label}/band"
[[ ! -e "${work}" ]] || exit 2
module purge
source /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS
export I_MPI_HYDRA_BOOTSTRAP=fork
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="${runtime}/source"
cd "${runtime}/source"
/public/home/iai806/.conda/envs/icu/bin/python "${root}/bundle/control.py" run-check \
 --bundle-sha256 "${BUNDLE_SHA256:?external bundle pin required}" --channel-id "${CHANNEL_ID}"
/public/home/iai806/.conda/envs/icu/bin/python -m vcneb.material_runner \
 --initial "${seed}/initial.vasp" --final "${seed}/final.vasp" \
 --initial-chain "${seed}/seed.traj" --n-images 9 --workdir "${work}" \
 --factory examples.hfo2_fixed_input_factory:make_clamped_seed_cached_factory \
 --parameters "${seed}/factory_parameters.json" \
 --command "mpirun -np 32 /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus" \
 --clamped-plane-reference "${seed}/substrate.vasp" --clamped-allow-tilt true --cell-scale 5.12968067458423 \
 --pressure-gpa 0 --cell-mode full --optimizer FIRE --fmax .10 --steps 60 --maxstep .02 --k .2 \
 --no-climb --mapping identity --mic --no-align-cells --cell-interpolation linear \
 --minimum-distance 1.6 --maximum-deformation .25 --maximum-cell-step .02 --candidate-step-retries 4 \
 --require-continuous-periodic-lift --summary-format brief
