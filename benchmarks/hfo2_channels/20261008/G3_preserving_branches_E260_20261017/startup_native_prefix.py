from pathlib import Path
from datetime import datetime,timezone,timedelta
import importlib.util,json,subprocess,hashlib
import numpy as np
root=Path('/public/home/iai806/abacus/agent-runs/20261017-varneb-G3-P-branches-E260-r2')
bundle=root/'bundle'
raw=subprocess.check_output(['sacct','-X','-n','-P','-j','28934048,28934049',
    '-o','JobIDRaw,State,ExitCode,ElapsedRaw,AllocCPUS,NodeList,Partition,NNodes,Start,End'],text=True)
paths=sorted(root.glob('shard_*/point_*/eval-*/result.json'))
failed=sorted(str(p) for p in root.glob('shard_*/point_*/eval-*/failure.json'))
points=[json.loads(p.read_bytes()) for p in sorted(root.glob('shard_*/point_*.json'))]
proof=dict(observed_CST=datetime.now(timezone(timedelta(hours=8))).isoformat(),
    raw_sacct=raw,complete_result_records_at_snapshot=len(paths),failure_records=failed,
    completed_point_rows=points,native_prefix=[],new_DFT_or_scheduler_writes=0,
    final_training_stability_or_prediction_certified=False)
if paths:
    spec=importlib.util.spec_from_file_location('E260_native_prefix_worker',bundle/'worker.py')
    worker=importlib.util.module_from_spec(spec);spec.loader.exec_module(worker)
    config=worker.legacy.pinned(bundle/'config.json','63b1aa0eea88e7cd7c2f1a78ce9f4f3c7f39e16a9caf04d3f6d8c6ef5f84ab80')
    reader,w,data,*_=worker.native_context(config,bundle)
    for shard in (0,1):
        shard_paths=[p for p in paths if p.relative_to(root).parts[0]==f'shard_{shard}']
        if not shard_paths:continue
        path=shard_paths[0];record=json.loads(path.read_bytes());x=np.asarray(record['coordinates'])
        native,source=reader.native_record(w,data,x,path,path.parent.parent)
        proof['native_prefix'].append(dict(shard=shard,record=native,source=source))
out=root/'first_native_prefix_2331.json'
with out.open('x',encoding='utf-8') as f:json.dump(proof,f,indent=2,allow_nan=False);f.write('\n')
print(json.dumps(dict(observed_CST=proof['observed_CST'],raw_sacct=raw,
    complete_result_records=len(paths),completed_points=len(points),failure_records=len(failed),
    audited_first_SCFs=[dict(shard=r['shard'],validation=r['record']['validation']) for r in proof['native_prefix']],
    report_sha256=hashlib.sha256(out.read_bytes()).hexdigest())))
