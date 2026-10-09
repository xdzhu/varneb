#!/bin/bash
# Same physical contract, twelve cached images, bounded fresh-state continuation.
#SBATCH -J hfo2-M-sample-refine
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=08:00:00
set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_JOB_NUM_NODES:-0}" == 1 && "${SLURM_NTASKS:-0}" == 1 && "${SLURM_CPUS_PER_TASK:-0}" == 32 ]] || exit 2
source_root=${SOURCE_ROOT:?set tested source}
seed_root=${SEED_ROOT:?set audited twelve-image seed}
workdir=${WORKDIR:?set fresh run directory}
[[ ! -e "${workdir}" ]] || exit 2
module purge
source /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS
export I_MPI_HYDRA_BOOTSTRAP=fork OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONPATH="${source_root}${PYTHONPATH:+:${PYTHONPATH}}"
cd "${source_root}"
/public/home/iai806/.conda/envs/icu/bin/python - "${seed_root}" <<'PY'
import json,sys,hashlib
from pathlib import Path
p=Path(sys.argv[1]); m=json.loads((p/'manifest.json').read_text())
assert m['purpose']=='G1_M_twelve_image_sampling_refinement' and m['n_total_images']==12
assert m['ordinary_target_eV_A']==.10 and m['step_cap']==20 and m['wall_cap_hours']==8
assert not m['climb'] and not m['physical_inputs_changed'] and m['all_twelve_initial_E_F_stress_reused']
assert all(hashlib.sha256((p/n).read_bytes()).hexdigest()==h for n,h in m['seed_file_sha256'].items())
PY
/public/home/iai806/.conda/envs/icu/bin/python -m vcneb.material_runner \
 --initial "${seed_root}/initial.vasp" --final "${seed_root}/final.vasp" --resume-snapshot "${seed_root}/seed.traj" \
 --n-images 12 --workdir "${workdir}" --factory examples.hfo2_fixed_input_factory:make_seed_cached_factory \
 --parameters "${seed_root}/factory_parameters.json" \
 --command "mpirun -np 32 /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus" \
 --pressure-gpa 0 --cell-mode full --optimizer FIRE --fmax .10 --steps 20 --maxstep .02 --k .2 --no-climb \
 --mapping identity --no-align-cells --cell-interpolation linear --minimum-distance 1.6 --maximum-deformation .25 \
 --maximum-cell-step .02 --candidate-step-retries 4 --require-continuous-periodic-lift --summary-format brief
