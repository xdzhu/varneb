set -euo pipefail
runtime=/public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-G2-E054-r1
prior=/public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-G2-E053-r1
[[ "$(sha256sum "$runtime/source.tar" | cut -d' ' -f1)" == 48c7c0e8550e3812d53d810c5e9d720f03693b0a7bc085856a663afc2e7712c2 ]] || exit 2
[[ ! -e "$runtime/observations" && ! -e "$runtime/audit_prepare.json" ]] || exit 2
for job in 28574708 28574709; do
 row=$(sacct -n -X -j "$job" --format=JobIDRaw,State,ExitCode -P | tr -d ' ')
 [[ "$row" == "$job|COMPLETED|0:0" ]] || { echo 'prior job not proved completed; stop' >&2; exit 2; }
done
export PYTHONPATH="$runtime/source" PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$runtime/source"
bash -n cluster/hf_hfo2_clamped_chain_resume_20261010.slurm
/public/home/iai806/.conda/envs/icu/bin/python - <<'PY'
import json
from pathlib import Path
import numpy as np
from ase.io import read
import examples.hfo2_fixed_input_factory as transport
from scripts.audit_hfo2_static_replica import sha256, audited_results
from scripts.export_hfo2_clamped_observation import export
from scripts.prepare_hfo2_clamped_resume import prepare
runtime=Path('/public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-G2-E054-r1')
prior=Path('/public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-G2-E053-r1')
def prohibit_external(*a,**k):
    raise RuntimeError('analysis/preflight prohibits external executable launch')
transport.subprocess.run=prohibit_external
reports=[]
for channel,job in [('PO_flip_T_pattern_preserving','28574708'),('PO_to_M','28574709')]:
    work=prior/channel/'band'
    summary=json.loads((work/'vcneb_summary.json').read_text())
    assert summary['status']=='max_steps_reached' and not summary['converged']
    assert not (work/'vcneb_failure.json').exists()
    observations=[]
    for step in (0,10):
        out=runtime/'observations'/channel/f'step_{step:04d}'
        r=export(work,step,out,job,production_script=prior/'source/cluster/hf_hfo2_clamped_chain_pilot_20261010.slurm')
        observations.append({k:r[k] for k in ('snapshot_step','replayed_fmax_eV_A','ordinary_residual_pass',
            'sampled_forward_barrier_meV_fu','sampled_reverse_barrier_meV_fu','evaluated_chain_sha256')})
    seed=runtime/'seeds'/channel
    m=prepare(runtime/'observations'/channel/'step_0010',seed)
    p=json.loads((seed/'factory_parameters.json').read_text())
    assert all(sha256(seed/n)==h for n,h in m['files_sha256'].items())
    images=read(seed/'seed.traj',index=':')
    factory=transport.make_clamped_resume_cached_factory(parameters=p,command=
       'mpirun -np 32 /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus')
    for i,a in enumerate(images):
        old=dict(a.calc.results)
        a.calc=factory(i,a,runtime/'cache_preflight'/channel/f'image_{i:04d}')
        for key,actual in [('energy',a.get_potential_energy()),('forces',a.get_forces()),('stress',a.get_stress())]:
            assert np.allclose(actual,old[key],atol=1e-12,rtol=0)
        assert a.calc.next_call==0
    calls=sorted(work.glob('image_*/scf_*/call_audit.json'))
    elapsed=0.
    for path in calls:
        audit=json.loads(path.read_text())
        assert {n:sha256(path.parent/n) for n in (*transport.CONTRACT,'STRU')}==audit['input_sha256']
        assert all(audit['input_sha256'][n]==h for n,h in transport.CONTRACT.items())
        assert sha256(path.parent/'OUT.ABACUS/running_scf.log')==audit['raw_log_sha256']
        raw=audited_results(path.parent)
        assert all(np.allclose(raw[k],audit['results'][k],atol=1e-12,rtol=0) for k in ('energy','forces','stress'))
        elapsed+=audit['elapsed_seconds']
    assert len(calls)==77
    reports.append({'channel':channel,'prior_job_id':job,'observations':observations,
       'fresh_prior_internal_SCFs':len(calls),'prior_SCF_transport_seconds_sum':elapsed,
       'prior_SCF_transport_core_hours_32':elapsed*32/3600,'all_original_six_physical_bytes_and_native_DSIZE32_checked':True,
       'current_frame_nine_exact_caches_checked':True,'external_launch_prohibited':True,'new_DFT_calls':0,
       'geometry_continuation_manifest_sha256':sha256(seed/'manifest.json'),
       'same_independent_chain':True,'FIRE_state_restored':False})
receipt={'status':'audited_clamped_geometry_continuations_prepared','reports':reports,
         'source_tree':'2c3aca9e49f75523c5b764d0e6554a5547a2d384','new_DFT_calls':0,
         'holdout_generated_or_read':False,'physical_inputs_changed':False,'ordinary_chains_converged':False}
(runtime/'audit_prepare.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))
PY
# End marker keeps PowerShell's final CRLF outside the here-document delimiter.
