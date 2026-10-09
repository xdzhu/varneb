set -euo pipefail
cd /public/home/iai806/abacus/agent-runs/20261009-varneb-clamped-PO-E046-r1/source
export PYTHONPATH="$PWD" OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
/public/home/iai806/.conda/envs/icu/bin/python - <<'PY'
import json
from pathlib import Path
import numpy as np
from ase.io import read
from examples.hfo2_fixed_input_factory import CONTRACT,read_fixed_hfo2_stru,same_ordered_geometry
from scripts.audit_hfo2_clamped_canary import replay_directories
from scripts.audit_hfo2_static_replica import sha256
from scripts.relax_clamped_ase_endpoint import load_seed
seed=Path("benchmarks/hfo2_channels/20261008/clamped_endpoint_seeds/strain_0000/PO_minus_T_preserving/endpoint_seed.json")
root=Path("/public/home/iai806/abacus/agent-runs/20261009-varneb-clamped-pair-E048-r1/PO_minus_T_preserving/endpoint")
assert sha256(Path("scripts/audit_hfo2_clamped_canary.py"))=="0a6700231b0086ff0fa6551f4b1b5b293c1a4d1f268590c534d7c60e11b67639"
assert sha256(seed)=="7d6131250346c6dea2565257f6a51b11acf5e8ca8fe29a08efd5241be383bf39"
atoms,boundary,m=load_seed(seed)
s=json.loads((root/"endpoint_relax_summary.json").read_text())
assert m["phase_label"]=="PO_minus_T_preserving" and m["strain"]==m["pressure_gpa"]==0
assert m["physical_contract_sha256"]==CONTRACT and s["physical_contract_sha256"]==CONTRACT
assert s["strain"]==s["external_pressure_gpa"]==0 and s["allow_tilt"] is True
assert s["steps_requested"]==20 and s["maxstep"]==.02 and s["fmax_target_eV_per_A"]==.03 and s["open_stress_target_kbar"]==2.
assert sha256(seed)==s["seed_manifest_sha256"]
calls=sorted((root/"calculator/image_0000").glob("scf_*"))
assert len(calls)==s["optimizer_steps"]+1 and 1<=len(calls)<=21
assert same_ordered_geometry(atoms,read_fixed_hfo2_stru(calls[0]/"STRU"))
rows=replay_directories(calls,boundary,full_physical_bytes=True)
assert same_ordered_geometry(read(root/"CONTCAR",format="vasp"),read_fixed_hfo2_stru(calls[-1]/"STRU"))
for k,v in (("potential_energy_eV",rows[-1]["energy_eV_cell"]),("max_atomic_force_eV_per_A",rows[-1]["max_atomic_force_eV_A"]),("open_traction_norm_kbar",rows[-1]["open_traction_norm_kbar"])):
    assert np.isclose(s[k],v,atol=2e-8,rtol=0)
physical=bool(s["max_atomic_force_eV_per_A"]<.03 and s["open_traction_norm_kbar"]<2.)
assert s["converged"] is physical and s["status"]==("completed" if physical else "step_limit")
report={"status":"raw_registered_uncached_endpoint_screen_passed_not_variant_or_G2_certification","job_id":"28456313","new_DFT_calls_for_analysis":0,"new_SCF_calls":len(rows),"BFGS_steps":s["optimizer_steps"],"endpoint_converged":s["converged"],"raw_full_physical_bytes_checked_here":True,"real_MPI_ranks":32,"seed_manifest_sha256":sha256(seed),"summary_sha256":sha256(root/"endpoint_relax_summary.json"),"replay_helper_sha256":sha256(Path("scripts/audit_hfo2_clamped_canary.py")),"SCF_seconds":sum(r["elapsed_seconds"] for r in rows),"SCF_core_hours":sum(r["elapsed_seconds"] for r in rows)*32/3600,"terminal_phase_symbols":[x["symbol"] for x in rows[-1]["structure_audit"]["symmetry_sweep"]],"rows":rows,"Hessian_variant_electronic_polarization_or_channel_barrier_certified":False,"original_source_or_runtime_overwritten":False}
print("E048_AUDIT_JSON="+json.dumps(report))
PY
sacct -j 28456313 --format=JobID,State,ExitCode,Start,End,Elapsed,AllocCPUS,NodeList -n -P
# End one completed endpoint zero-DFT raw audit.
