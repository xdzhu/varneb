"""Native read-only E161 prefix audit plus exact cache-only hole attribution."""
from pathlib import Path
from types import SimpleNamespace
import argparse
import hashlib
import importlib.util
import json
import sys

sys.dont_write_bytecode=True


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--driver-sha256',required=True)
    p.add_argument('--diagnostics-sha256',required=True)
    args=p.parse_args()
    if args.output.exists():raise FileExistsError('fresh external report required')
    output_root=Path('/public/home/iai806/abacus/agent-runs/20261015-varneb-G3-bound-audit-E168-r1')
    if args.output.resolve()!=output_root/'preserving_native_prefix_and_holes.json':
        raise ValueError('fresh dedicated E168 report path required')
    sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
    local=Path(__file__).resolve().parent
    if sha(Path(__file__))!=args.driver_sha256 or sha(local/'diagnostics.py')!=args.diagnostics_sha256:
        raise ValueError('diagnostic source bytes differ')
    old=Path('/public/home/iai806/abacus/agent-runs/20261015-varneb-G3-face-audit-E163-r1/audit_hfo2_G3_conditional_faces.py')
    if sha(old)!='ab041335a93b8eaa98fa5cf8a812db5126931105cf972f586e4d15f8d06742da':
        raise ValueError('existing tested auditor bytes differ')
    import diagnostics
    spec=importlib.util.spec_from_file_location('original_E163_native_auditor',old)
    audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
    original=audit.replay_node;holes=[]
    proxy=SimpleNamespace(worker=audit.worker,replay_node=original)
    def observed(data,row,records,start):
        if row['point'] is not None:return original(data,row,records,start)
        detail=diagnostics.observe_replay(proxy,data,row,records,start)
        holes.append(detail)
        return detail['replay']
    try:
        audit.replay_node=observed
        report=audit.audit('preserving','prefix')
    finally:
        audit.replay_node=original
    expected={r['id'] for r in report['candidate_rows'] if r['candidate'] is None and r['completed_SCFs']>0}
    if {r['point_id'] for r in holes}!=expected or 13 not in expected:
        raise ValueError('all completed nonconverged nodes including original point13 must be diagnosed')
    report.update(cache_only_hole_diagnostics=holes,diagnostic_source_sha256=args.diagnostics_sha256,
        native_driver_sha256=args.driver_sha256,optimizer_or_domain_or_input_changes=0,
        full_goal_active_incomplete=True)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    audit.study.exclusive_json(args.output,report)
    print(json.dumps(dict(status=report['status'],completed_nodes=report['completed_nodes_audited'],
        raw_SCFs_reparsed=report['new_SCFs_reparsed'],holes=[dict(point_id=r['point_id'],
            guard_counts=r['guard_counts'],last_trial=r['completed_trial_rows'][-1]) for r in holes],
        new_DFT=0,report_sha256=sha(args.output))))


if __name__=='__main__':main()
