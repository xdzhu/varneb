#!/bin/bash
# Bounded ordinary continuation after a terminal, exact-SCF health segment.
#SBATCH -J hfo2-g1-continue
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_JOB_NUM_NODES:-0}" == 1 \
   && "${SLURM_NTASKS:-0}" == 1 && "${SLURM_CPUS_PER_TASK:-0}" == 32 ]] || exit 2
source_root=${SOURCE_ROOT:?set immutable production archive}
seed_root=${SEED_ROOT:?set terminal-run exact-cache restart}
workdir=${WORKDIR:?set fresh continuation workdir}
[[ ! -e "${workdir}" ]] || { echo "refusing existing run" >&2; exit 2; }
module purge
source /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS
export I_MPI_HYDRA_BOOTSTRAP=fork
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONPATH="${source_root}${PYTHONPATH:+:${PYTHONPATH}}"
cd "${source_root}"
n_images=$(/public/home/iai806/.conda/envs/icu/bin/python - "${seed_root}" <<'PY'
import hashlib,json,sys
from pathlib import Path
p=Path(sys.argv[1]); m=json.loads((p/'manifest.json').read_text())
assert m['purpose']=='ordinary_observation_restart'
assert m['n_total_images'] in (9,10) and m['fmax_eV_A']==0.10
assert not m['climb'] and not m['physical_inputs_changed']
assert m['all_seed_SCFs_reused'] and not m['old_optimizer_state_reused']
assert m['next_optimizer_step_cap']==80 and m['wall_cap_hours']==24
assert all(hashlib.sha256((p/n).read_bytes()).hexdigest()==h for n,h in m['seed_file_sha256'].items())
print(m['n_total_images'])
PY
)
/public/home/iai806/.conda/envs/icu/bin/python -m vcneb.material_runner \
  --initial "${seed_root}/initial.vasp" --final "${seed_root}/final.vasp" \
  --resume-snapshot "${seed_root}/seed.traj" --n-images "${n_images}" --workdir "${workdir}" \
  --factory examples.hfo2_fixed_input_factory:make_seed_cached_factory \
  --parameters "${seed_root}/factory_parameters.json" \
  --command "mpirun -np 32 /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus" \
  --pressure-gpa 0 --cell-mode full --optimizer FIRE --fmax 0.10 --steps 80 \
  --maxstep 0.02 --k 0.2 --no-climb --mapping identity --mic \
  --no-align-cells --cell-interpolation linear --minimum-distance 1.6 \
  --maximum-deformation 0.25 --maximum-cell-step 0.02 --candidate-step-retries 4 \
  --require-continuous-periodic-lift --summary-format brief
