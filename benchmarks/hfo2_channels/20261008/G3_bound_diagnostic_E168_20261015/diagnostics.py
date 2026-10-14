"""Cache-only guard attribution; never alter decisions or request new E/g."""
from __future__ import annotations
from collections import Counter
import hashlib
import numpy as np


def guard_reason(chart,basis,coordinates,settings):
    """Follow the original guard's short-circuit order, not an alternate guard."""
    x=np.asarray(coordinates,dtype=float)
    if x.shape!=(36,):raise ValueError('original 36D internal coordinates required')
    if not np.isfinite(x).all():
        return dict(safe=False,reason='nonfinite',internal_norm_A=None,minimum_distance_A=None)
    norm=float(np.linalg.norm(x))
    if norm>settings['maximum_internal_displacement_A']:
        return dict(safe=False,reason='original_internal_radius',internal_norm_A=norm,minimum_distance_A=None)
    try:
        atoms=chart.to_atoms(basis@x)
        distances=atoms.get_all_distances(mic=True)
        np.fill_diagonal(distances,np.inf)
        minimum=float(np.min(distances))
        safe=bool(minimum>=settings['minimum_distance_A'])
        return dict(safe=safe,reason='accepted' if safe else 'minimum_distance',
                    internal_norm_A=norm,minimum_distance_A=minimum)
    except (ValueError,np.linalg.LinAlgError):
        return dict(safe=False,reason='chart_or_cell_invalid',internal_norm_A=norm,minimum_distance_A=None)


def trial_rows(data,records,radius,target):
    """Every completed trial, not a claim about accepted optimizer iterates."""
    basis,axes=data['basis'],data['axes'];rows=[]
    for i,record in enumerate(records):
        x=basis.T@np.asarray(record['coordinates']);g=basis.T@np.asarray(record['gradient'])
        q=axes.T@x;r=x-axes@q;gr=g-axes@(axes.T@g)
        if x.shape!=(36,) or not np.isfinite(x).all() or not np.isfinite(gr).all():
            raise ValueError('finite original trial pullback required')
        rnorm=float(np.linalg.norm(r));norm=float(np.linalg.norm(x));product=float(r@gr)
        rows.append(dict(trial_index=i,internal_norm_A=norm,radius_margin_A=radius-norm,
            q_A=q.tolist(),released_norm_A=rnorm,released_gradient_norm_eV_A=float(np.linalg.norm(gr)),
            radius_squared_steepest_descent_derivative_eV=-2*product,
            steepest_descent_points_outward=bool(product<0),
            energy_change_from_original_center_meV_cell=1000*(record['enthalpy_eV']-data['center_energy']),
            below_gradient_target=bool(np.linalg.norm(gr)<=target),
            accepted_iterate_or_final_candidate_inferred=False))
    if not rows:raise ValueError('at least one complete native trial required')
    return rows


def observe_replay(audit,data,row,records,start):
    """Original replay checks order, counters and stop reason; restore hook even on error."""
    canary=audit.worker.canary;original=canary.geometry_guard;events=[]
    def observe(chart,basis,x):
        decision=original(chart,basis,x)
        detail=guard_reason(chart,basis,x,canary.SETTINGS)
        if detail['safe']!=decision:raise ValueError('guard classifier differs from original decision')
        events.append(dict(index=len(events),coordinate_sha256=hashlib.sha256(np.asarray(x).tobytes()).hexdigest(),**detail))
        return decision
    try:
        canary.geometry_guard=observe
        replay=audit.replay_node(data,row,records,start)
    finally:
        canary.geometry_guard=original
    counts=dict(Counter(e['reason'] for e in events))
    if (len(events)!=row['geometry_screens']['checks']
            or sum(not e['safe'] for e in events)!=row['geometry_screens']['rejected']):
        raise ValueError('guard event counters differ')
    trials=trial_rows(data,records,canary.SETTINGS['maximum_internal_displacement_A'],
                      canary.SETTINGS['gradient_tolerance_eV_A'])
    return dict(point_id=row['id'],q_A=row['q_A'],status=row['status'],stop_reason=row.get('reason'),
        replay=replay,guard_counts=counts,guard_events=events,completed_trial_rows=trials,
        completed_trials=len(trials),reported_final_candidate=row['point'],
        last_trial_not_promoted=True,new_DFT_calls=0,optimizer_or_guard_decisions_changed=False,
        unconstrained_minimum_unreachable_or_boundary_KKT_or_branch_stability_certified=False)
