#!/bin/bash
# Same registered flip chain: latest complete E067 step20 (cumulative50).
# Only the finite update allowance extends to40; all physical inputs unchanged.
# Geometry continuation resets FIRE history, not an optimizer checkpoint.
#SBATCH -J hfo2-G2-flip-E069
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=12:00:00
set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_JOB_NUM_NODES:-0}" == 1 \
   && "${SLURM_NTASKS:-0}" == 1 && "${SLURM_CPUS_PER_TASK:-0}" == 32 ]] || exit 2
source_root=/public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-G2-E054-r1/source-fixed
seed_root=${SEED_ROOT:?set audited latest step20 flip continuation}
seed_manifest_sha=${SEED_MANIFEST_SHA256:?set exact audited manifest digest}
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
 'import hashlib,json,pathlib,sys
p=pathlib.Path(sys.argv[1]); digest=sys.argv[2]
if len(digest)!=64 or hashlib.sha256(p.read_bytes()).hexdigest()!=digest: raise ValueError("seed manifest changed")
m=json.loads(p.read_text())
required={"status":"audited_G2_geometry_continuation_prepared","source_job_id":"28722810","source_step":20,"current_frame_exact_caches":9,"new_independent_chains":0,"new_DFT_calls":0,"physical_inputs_changed":False,"climb":False,"holdout_generated_or_read":False,"optimizer_state_restored":False,"pressure_GPa":0,"ordinary_fmax_eV_A":.10,"n_total_images":9,"n_fixed_endpoints":2,"n_active_images":7,"cell_scale_A":5.12968067458423}
if any(m.get(k)!=v for k,v in required.items()): raise ValueError("registered same-chain seed contract differs")
names={"initial.vasp","final.vasp","seed.traj","substrate.vasp","factory_parameters.json"}
if set(m["files_sha256"])!=names: raise ValueError("seed file set changed")
if any(hashlib.sha256((p.parent/n).read_bytes()).hexdigest()!=h for n,h in m["files_sha256"].items()): raise ValueError("seed file changed")
print(m["cell_scale_A"])' "${seed_root}/manifest.json" "${seed_manifest_sha}")
/public/home/iai806/.conda/envs/icu/bin/python -m vcneb.material_runner \
 --initial "${seed_root}/initial.vasp" --final "${seed_root}/final.vasp" \
 --initial-chain "${seed_root}/seed.traj" --n-images 9 --workdir "${workdir}" \
 --factory examples.hfo2_fixed_input_factory:make_clamped_resume_cached_factory \
 --parameters "${seed_root}/factory_parameters.json" \
 --command "mpirun -np 32 /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus" \
 --clamped-plane-reference "${seed_root}/substrate.vasp" --clamped-allow-tilt true --cell-scale "${cell_scale}" \
 --pressure-gpa 0 --cell-mode full --optimizer FIRE --fmax .10 --steps 40 --maxstep .02 --k .2 \
 --no-climb --mapping identity --mic --no-align-cells --cell-interpolation linear \
 --minimum-distance 1.6 --maximum-deformation .25 --maximum-cell-step .02 --candidate-step-retries 4 \
 --require-continuous-periodic-lift --summary-format brief
