"""Arithmetic from the pinned native replay, not a new optimization or DFT."""
from pathlib import Path
import argparse
from collections import Counter
import hashlib
import json
import math

CASE=Path(__file__).resolve().parent
NATIVE_SHA='7c9581eb0c102f6110415ca590b7f045b4875e3d3498a7af5019d0cb8a14f57b'


def radial_split(released_radius,gradient_norm,squared_radius_descent_derivative):
    values=(released_radius,gradient_norm,squared_radius_descent_derivative)
    if not all(math.isfinite(v) for v in values) or released_radius<=0 or gradient_norm<0:
        raise ValueError('finite positive released radius and nonnegative gradient norm required')
    radial=-squared_radius_descent_derivative/(2*released_radius)
    square=gradient_norm**2-radial**2
    if square < -1e-12*max(gradient_norm**2,1e-30):
        raise ValueError('radial component exceeds total gradient; not clipped')
    # Only possible negative roundoff in this scalar norm identity, not curvature repair.
    return dict(radial_gradient_component_eV_A=radial,
                constant_radius_tangent_gradient_norm_eV_A=math.sqrt(max(square,0.)))


def analyze():
    path=CASE/'preserving_native_prefix_and_holes.json';raw=path.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=NATIVE_SHA:raise ValueError('native report bytes differ')
    source=json.loads(raw)
    if (source['job_id']!='28862441' or source['channel']!='preserving'
            or source['new_DFT_or_scheduler_writes']!=0
            or source['runtime_manifest_sha256']!='d07bb011c9e313afe23ca4e78207356ff98adfc43ca231bdd047e6d5d0f1be20'
            or source['auditor_sha256']!='ab041335a93b8eaa98fa5cf8a812db5126931105cf972f586e4d15f8d06742da'):
        raise ValueError('original runtime/auditor/job identity differs')
    holes=[]
    for hole in source['cache_only_hole_diagnostics']:
        if (hole['reported_final_candidate'] is not None or hole['optimizer_or_guard_decisions_changed']
                or hole['replay']!=dict(original_order_cache_only_replay_passed=True,new_DFT=0)):
            raise ValueError('original hole/replay identity differs')
        counts=dict(Counter(e['reason'] for e in hole['guard_events']))
        if counts!=hole['guard_counts']:raise ValueError('guard event counts differ')
        rows=[dict(row,**radial_split(row['released_norm_A'],row['released_gradient_norm_eV_A'],
            row['radius_squared_steepest_descent_derivative_eV'])) for row in hole['completed_trial_rows']]
        if len(rows)!=12 or any(r['accepted_iterate_or_final_candidate_inferred'] for r in rows):
            raise ValueError('original twelve trial records required; no candidate promotion')
        holes.append(dict(point_id=hole['point_id'],q_A=hole['q_A'],guard_counts=counts,
            stop_reason=hole['stop_reason'],trial_rows=rows,
            minimum_observed_trial_gradient_norm_eV_A=min(r['released_gradient_norm_eV_A'] for r in rows),
            last_completed_trial=rows[-1],original_radius_A=.5,original_gradient_target_eV_A=.03,
            no_completed_trial_below_target=not any(r['below_gradient_target'] for r in rows),
            boundary_KKT_or_unreachable_unconstrained_minimum_or_unique_cause_certified=False))
    return dict(status='same_order_native_hole_replay_with_scalar_radial_diagnostic',
        native_report_sha256=NATIVE_SHA,original_prefix_nodes=source['completed_nodes_audited'],
        original_prefix_candidates=source['gradient_candidates'],raw_SCFs_reparsed_by_native_driver=source['new_SCFs_reparsed'],
        diagnostic_source_sha256=source['diagnostic_source_sha256'],native_driver_sha256=source['native_driver_sha256'],
        holes=holes,radial_split_is_fixed_metric_geometry_not_phonon_or_barrier_energy=True,
        rejected_outside_radius_proposals_not_distance_tested=True,
        no_complete_PES_TS_or_independent_material_prediction_certified=True,
        new_DFT_or_scheduler_writes=0,new_raw_SCF_reparse_in_this_arithmetic=False,
        original_input_domain_force_target_and_runtime_changed=False,
        physical_reserved_geometry_generated_or_read=False)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if args.output.exists():raise FileExistsError('fresh arithmetic output required')
    result=analyze();args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x',encoding='utf-8',newline='\n') as f:
        json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(dict(status=result['status'],holes=len(result['holes']),new_DFT=0)))


if __name__=='__main__':main()
