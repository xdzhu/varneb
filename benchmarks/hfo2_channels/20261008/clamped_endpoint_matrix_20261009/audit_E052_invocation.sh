set -euo pipefail
runtime=/public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-p0100-minus-E052-r1
source_root=/public/home/iai806/abacus/agent-runs/20261009-varneb-clamped-PO-E046-r1/source
cd "$source_root"
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
runtime=Path('/public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-p0100-minus-E052-r1')
for phase,job in [('PO_minus_T_preserving','28568571'),('PO_minus_T_reversing','28568573')]:
    seed=Path('benchmarks/hfo2_channels/20261008/clamped_endpoint_seeds/strain_p0100')/phase/'endpoint_seed.json'
    root=runtime/phase/'endpoint'
    output=runtime/phase/'raw_audit_HF.json'
    assert not output.exists(), 'fresh audit artifact required'
    initial,boundary,m=load_seed(seed)
    s=json.loads((root/'endpoint_relax_summary.json').read_text())
    assert m['phase_label']==phase and m['strain']==s['strain']==.01 and m['pressure_gpa']==s['external_pressure_gpa']==0
    assert m['physical_contract_sha256']==s['physical_contract_sha256']==CONTRACT
    assert s['allow_tilt'] is True and s['steps_requested']==20 and s['maxstep']==.02
    assert s['fmax_target_eV_per_A']==.03 and s['open_stress_target_kbar']==2.
    assert s['seed_manifest_sha256']==sha256(seed)
    calls=sorted((root/'calculator/image_0000').glob('scf_*'))
    assert 1<=len(calls)<=21 and len(calls)==s['optimizer_steps']+1
    assert same_ordered_geometry(initial,read_fixed_hfo2_stru(calls[0]/'STRU'))
    rows=replay_directories(calls,boundary,full_physical_bytes=True)
    assert same_ordered_geometry(read(root/'CONTCAR',format='vasp'),read_fixed_hfo2_stru(calls[-1]/'STRU'))
    for k,v in [('potential_energy_eV',rows[-1]['energy_eV_cell']),('max_atomic_force_eV_per_A',rows[-1]['max_atomic_force_eV_A']),('open_traction_norm_kbar',rows[-1]['open_traction_norm_kbar'])]:
        assert np.isclose(s[k],v,atol=2e-8,rtol=0)
    physical=bool(s['max_atomic_force_eV_per_A']<.03 and s['open_traction_norm_kbar']<2.)
    assert physical and s['converged'] is physical and s['status']=='completed'
    assert [p['symbol'] for p in rows[-1]['structure_audit']['symmetry_sweep']]==['Pca2_1']*3
    report={'status':'raw_registered_uncached_endpoint_screen_passed_not_variant_or_G2_certification','job_id':job,'new_DFT_calls_for_analysis':0,'new_SCF_calls':len(rows),'BFGS_steps':s['optimizer_steps'],'endpoint_converged':s['converged'],'raw_full_physical_bytes_checked_here':True,'real_MPI_ranks':32,'seed_manifest_sha256':sha256(seed),'summary_sha256':sha256(root/'endpoint_relax_summary.json'),'SCF_seconds':sum(r['elapsed_seconds'] for r in rows),'SCF_core_hours':sum(r['elapsed_seconds'] for r in rows)*32/3600,'terminal_phase_symbols':[p['symbol'] for p in rows[-1]['structure_audit']['symmetry_sweep']],'rows':rows,'Hessian_variant_electronic_polarization_or_channel_barrier_certified':False,'original_source_or_runtime_overwritten':False}
    output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='rows'}))
PY
for phase in PO_minus_T_preserving PO_minus_T_reversing; do
 cd "$runtime/$phase"
 [[ ! -e completed_observables.tar.gz ]] || exit 2
 find endpoint -type f \( -name INPUT -o -name KPT -o -name STRU -o -name 'abacus.out' -o -name 'abacus.err' -o -name 'input_sha256.json' -o -name 'call_audit.json' -o -name 'running_scf.log' -o -name 'relax.log' -o -name 'relax.traj' -o -name 'CONTCAR*' -o -name 'endpoint_relax*.json' -o -name 'structure.start.vasp' \) -print0 | tar --null -czf completed_observables.tar.gz -T -
 sha256sum raw_audit_HF.json completed_observables.tar.gz
done
sacct -j 28568571,28568573 --format=JobID,State,ExitCode,Start,End,Elapsed,AllocCPUS,NodeList -n -P
# End zero-DFT audit and observable-only export, no licensed basis or charges.
