#!/bin/bash
# Finite continuation of the audited +1% reversing child40 cap.
# Latest complete cap audit and all nine exact caches are required; no dispatch.
# The E090 driver/settings are unchanged; only registered ancestry differs.
#SBATCH -J hfo2-G2-continue-E146
#SBATCH -p hfacnormal01
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=12:00:00
set -euo pipefail
[[ "${RUN_DFT:-0}" == 1 && "${SLURM_JOB_NUM_NODES:-0}" == 1 \
   && "${SLURM_NTASKS:-0}" == 1 && "${SLURM_CPUS_PER_TASK:-0}" == 32 ]] || exit 2
source_root=${SOURCE_ROOT:?set unchanged original parent runtime}
seed_root=${SEED_ROOT:?set actual audited latest complete nine-cache seed}
workdir=${WORKDIR:?set fresh optimization namespace}
parent_job=${PARENT_JOB_ID:?set registered terminal parent}
parent_audit=${PARENT_AUDIT_PATH:?set actual reviewed terminal cap receipt}
parent_audit_sha=${PARENT_AUDIT_SHA256:?set exact reviewed receipt digest}
seed_sha=${SEED_MANIFEST_SHA256:?set exact reviewed latest-seed digest}
[[ ! -e "${workdir}" ]] || exit 2
cell_scale=$(/public/home/iai806/.conda/envs/icu/bin/python -c \
 'import hashlib,json,pathlib,sys
seed,seed_digest,audit,audit_digest,parent,source=sys.argv[1:]
runs="/public/home/iai806/abacus/agent-runs/"
registered={"28764704":(40,0,338,runs+"20261010-varneb-nested-controls-E055-r1/source",runs+"20261012-varneb-clamped-flip-continue-E103-r1/28709791/PO_flip_T_pattern_reversing/band")}
if parent not in registered: raise ValueError("only the reviewed reversing child40 parent is admitted")
step,initial,code_count,old_source,old_work=registered[parent]
if source!=old_source: raise ValueError("unchanged original parent runtime required")
def pinned(path,digest):
 p=pathlib.Path(path)
 if len(digest)!=64 or hashlib.sha256(p.read_bytes()).hexdigest()!=digest: raise ValueError("reviewed input digest differs")
 return json.loads(p.read_text())
r=pinned(audit,audit_digest); m=pinned(seed,seed_digest)
required={"status":"audited_G2_geometry_continuation_prepared","source_job_id":parent,"source_step":step,"current_frame_exact_caches":9,"new_independent_chains":0,"new_DFT_calls":0,"physical_inputs_changed":False,"climb":False,"holdout_generated_or_read":False,"optimizer_state_restored":False,"pressure_GPa":0,"ordinary_fmax_eV_A":.10,"n_total_images":9,"n_fixed_endpoints":2,"n_active_images":7,"cell_scale_A":5.12968067458423}
if any(type(m.get(k)) is not type(v) or m.get(k)!=v for k,v in required.items()): raise ValueError("latest same-chain seed contract differs")
if r.get("status")!="audited_step_cap_not_converged_same_chain_seed" or r.get("job_id")!=parent or r.get("workdir")!=old_work or r.get("terminal_step")!=step or r.get("ordinary_converged") is not False or r.get("scheduler_row")!=parent+"|COMPLETED|0:0|32|hfacnormal01|iai806": raise ValueError("actual completed unconverged step-cap parent required")
if r.get("fresh_interior_SCFs")!=7*step+initial or r.get("initial_fresh_interior_SCFs")!=initial or r.get("current_frame_exact_caches_checked")!=9 or r.get("seed_manifest_sha256")!=seed_digest or r.get("all_calls_original_six_physical_bytes_native_DSIZE32_and_raw_EFS_checked") is not True: raise ValueError("all native SCFs and nine latest caches must be audited")
if any(r.get(k) is not False for k in ("physical_inputs_or_running_source_changed","holdout_generated_or_read","FIRE_state_restored")) or any(r.get(k)!=0 for k in ("new_DFT_calls","scheduler_mutations","new_independent_chains")): raise ValueError("audit scope differs")
obs=pathlib.Path(audit).parent/"observations"/("step_%04d"%step)/"observation.json"
if hashlib.sha256(obs.read_bytes()).hexdigest()!=m.get("source_observation_sha256") or json.loads(obs.read_text())!=r.get("terminal_observation"): raise ValueError("latest observation/receipt identity differs")
names={"initial.vasp","final.vasp","seed.traj","substrate.vasp","factory_parameters.json"}
if set(m.get("files_sha256",{}))!=names or any(hashlib.sha256((pathlib.Path(seed).parent/n).read_bytes()).hexdigest()!=h for n,h in m["files_sha256"].items()): raise ValueError("one of five exact seed files changed")
code=r.get("runtime_code_sha256",{})
if len(code)!=code_count or "vcneb/material_runner.py" not in code or any(hashlib.sha256((pathlib.Path(source)/n).read_bytes()).hexdigest()!=h for n,h in code.items()): raise ValueError("immutable parent runtime changed")
print(m["cell_scale_A"])' "${seed_root}/manifest.json" "${seed_sha}" "${parent_audit}" "${parent_audit_sha}" "${parent_job}" "${source_root}")
module purge
source /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/scripts/env.sh
unset I_MPI_PMI_LIBRARY I_MPI_HYDRA_BOOTSTRAP_EXEC_EXTRA_ARGS
export I_MPI_HYDRA_BOOTSTRAP=fork
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONPATH="${source_root}${PYTHONPATH:+:${PYTHONPATH}}"
cd "${source_root}"
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
