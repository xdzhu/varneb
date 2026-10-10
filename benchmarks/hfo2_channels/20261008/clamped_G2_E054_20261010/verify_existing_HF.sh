set -euo pipefail
runtime=/public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-G2-E054-r1
[[ ! -e "$runtime/verification_after_transport_error.json" ]] || exit 2
export PYTHONPATH="$runtime/source" PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
/public/home/iai806/.conda/envs/icu/bin/python - <<'PY'
import json
from pathlib import Path
import numpy as np
from ase.io import read
import examples.hfo2_fixed_input_factory as transport
from scripts.audit_hfo2_static_replica import audited_results,sha256
runtime=Path('/public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-G2-E054-r1')
def prohibit(*a,**k):
    raise RuntimeError('verification must not launch any external executable')
transport.subprocess.run=prohibit
receipt=json.loads((runtime/'audit_prepare.json').read_text())
assert receipt['status']=='audited_clamped_geometry_continuations_prepared' and receipt['new_DFT_calls']==0
assert len(receipt['reports'])==2 and receipt['source_tree']=='2c3aca9e49f75523c5b764d0e6554a5547a2d384'
checks=[]
for r in receipt['reports']:
    assert r['fresh_prior_internal_SCFs']==77 and r['all_original_six_physical_bytes_and_native_DSIZE32_checked']
    assert r['current_frame_nine_exact_caches_checked'] and r['external_launch_prohibited']
    seed=runtime/'seeds'/r['channel']
    assert sha256(seed/'manifest.json')==r['geometry_continuation_manifest_sha256']
    m=json.loads((seed/'manifest.json').read_text())
    assert all(sha256(seed/n)==h for n,h in m['files_sha256'].items())
    p=json.loads((seed/'factory_parameters.json').read_text())
    images=read(seed/'seed.traj',index=':')
    assert len(images)==len(p['seed_cache_records'])==9
    for i,(a,c) in enumerate(zip(images,p['seed_cache_records'])):
        source=Path(c['directory'])
        assert {n:sha256(source/n) for n in (*transport.CONTRACT,'STRU')}==c['input_sha256']
        assert all(c['input_sha256'][n]==h for n,h in transport.CONTRACT.items())
        assert sha256(source/'OUT.ABACUS/running_scf.log')==c['raw_log_sha256']
        assert transport.same_ordered_geometry(a,transport.read_fixed_hfo2_stru(source/'STRU'))
        raw=audited_results(source)
        assert all(np.allclose(raw[k],a.calc.results[k],atol=1e-12,rtol=0) for k in ('energy','forces','stress'))
        cached=json.loads((runtime/'cache_preflight'/r['channel']/f'image_{i:04d}/seed_cache_audit.json').read_text())
        assert cached['policy']=='identical_ordered_clamped_resume_hash_pinned' and cached['new_DFT_calls']==0
        assert cached['raw_source']==c['directory'] and cached['input_sha256']==c['input_sha256']
        assert cached['raw_log_sha256']==c['raw_log_sha256']
    checks.append({'channel':r['channel'],'exact_geometry_native_raw_EFS_and_pins_checked':9})
result={'status':'passed','audit_prepare_sha256':sha256(runtime/'audit_prepare.json'),
        'verified_complete_cached_frames':checks,'external_launch_prohibited':True,'new_DFT_calls':0,
        'original_harness_exit':1,'original_harness_error':'Final PowerShell CRLF became PY token after completed audit receipt',
        'repair':'Add final shell comment after here-doc; validate existing complete artifacts without re-export or cache rewrite',
        'original_source_raw_data_or_DFT_parameters_modified':False}
(runtime/'verification_after_transport_error.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
PY
# Explicit footer preserves the LF here-document delimiter through PowerShell.
