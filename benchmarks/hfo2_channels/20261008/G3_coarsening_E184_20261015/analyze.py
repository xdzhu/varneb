"""Post-hoc energy-only coarsening, not independent material validation.

Standard bilinear/tensor-quadratic interpolation and a refitted quadratic
provide comparators, not innovations. Never fill a missing nonzero-weight
node or label an optimizer trial as a converged candidate.
"""
from pathlib import Path
import argparse,hashlib,itertools,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT))
from scripts import hfo2_G3_conditional_faces as faces

SOURCE=ROOT/'benchmarks/hfo2_channels/20261008/G3_shared_Q_controls_E180_20261015/analysis.json'
PIN='39e31d790592e7a2f640a8428f00efecc28fd14775b7fb4b2b699b7da54c6402'
PROOFS={'T':('G3_feasible_execution_E171_20261015/T_complete_native_audit.json',
    'b8e7836722e3c7cadfaed9771eeb90e1d53097a03fc5d19272109eca0a47c166'),
    'preserving':('G3_preserving_terminal_E179_20261015/preserving_complete_native_audit.json',
    '982d2657a7a2c536509d5a8c460c1ccc0081028b367de398847720fe46937d5a')}
MODELS=['bilinear','tensor_quadratic','refitted_total_quadratic','original_H_0.01','original_H_0.02']


def pinned(path,digest):
    raw=path.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=digest:raise ValueError('original source bytes differ')
    return json.loads(raw)


def one_axis(index,kind):
    if index not in (-2,-1,0,1,2):raise ValueError('original registered grid indices only')
    x=index/2.
    if kind=='tensor_quadratic':pairs=zip((-2,0,2),(.5*x*(x-1),1-x*x,.5*x*(x+1)))
    elif kind=='bilinear':
        pairs=zip((-2,0),(-x,1+x)) if index<=0 else zip((0,2),(1-x,x))
    else:raise ValueError('declared elementary interpolation only')
    return {i:float(v) for i,v in pairs if v!=0.}


def interpolate(index,labels,kind):
    wx,wy=(one_axis(i,kind) for i in index)
    weights={(i,j):u*v for (i,u),(j,v) in itertools.product(wx.items(),wy.items())}
    missing=[list(ij) for ij in weights if ij not in labels]
    if missing:return dict(value=None,missing_required_nodes=missing,required_weights=[dict(index=list(k),weight=v) for k,v in weights.items()])
    value=sum(w*labels[ij] for ij,w in weights.items())
    return dict(value=float(value),missing_required_nodes=[],required_weights=[dict(index=list(k),weight=v) for k,v in weights.items()])


def features(index):
    x,y=np.asarray(index,dtype=float)/2.
    return np.array([1.,x,y,x*x,x*y,y*y])


def fit_quadratic(labels):
    keys=sorted(labels);design=np.array([features(k) for k in keys]);target=np.array([labels[k] for k in keys])
    rank=int(np.linalg.matrix_rank(design))
    if rank<6:return None,rank  # no unidentified coefficient or pseudoinverse fallback
    beta,_,returned_rank,_=np.linalg.lstsq(design,target,rcond=None)
    if returned_rank!=6:raise ValueError('quadratic fit rank changed')
    return beta,rank


def statistics(points,ids,model):
    selected=[p for p in points if p['id'] in ids]
    if any(p['models'][model]['value'] is None or p['observed'] is None for p in selected):
        raise ValueError('statistics require supported predictions and actual candidate labels')
    raw=[p['models'][model]['raw_error_meV_cell'] for p in selected]
    shape=[p['models'][model]['shape_error_meV_cell'] for p in selected]
    return dict(ids=[p['id'] for p in selected],count=len(selected),
        max_absolute_raw_error_meV_cell=None if not selected else max(map(abs,raw)),
        max_absolute_seen_origin_shape_error_meV_cell=None if not selected else max(map(abs,shape)),
        seen_origin_shape_RMSE_meV_cell=None if not selected else float(np.linalg.norm(shape)/np.sqrt(len(shape))))


def analyze():
    source=pinned(SOURCE,PIN);specs=faces.ordered_points();result=[]
    coarse={s['id'] for s in specs if all(i%2==0 for i in s['index'])}
    if len(coarse)!=9:raise ValueError('fixed geometric parity split required')
    for channel in ('T','preserving'):
        rawrows=[next(r for r in source['rows'] if r['channel']==channel and r['reference']=='T_all_optical'
                     and r['step_A']==h) for h in (.01,.02)]
        if any([p['q_A'] for p in r['points']]!=[s['q_A'] for s in specs] for r in rawrows):
            raise ValueError('original shared-observable ordering differs')
        observed={p['id']:p['DFT_energy_change_eV_cell'] for p in rawrows[0]['points']
                  if p['DFT_energy_change_eV_cell'] is not None}
        labels={tuple(s['index']):observed[s['id']] for s in specs if s['id'] in coarse and s['id'] in observed}
        beta,rank=fit_quadratic(labels)
        origins={kind:interpolate([0,0],labels,kind)['value'] for kind in ('bilinear','tensor_quadratic')}
        origins['refitted_total_quadratic']=None if beta is None else float(beta[0])
        origins.update({f'original_H_{h:.2f}':r['points'][0]['model_energy_change_eV_cell'] for h,r in zip((.01,.02),rawrows)})
        points=[]
        for s in specs:
            if s['id'] in coarse:continue
            predictions={kind:interpolate(s['index'],labels,kind) for kind in ('bilinear','tensor_quadratic')}
            predictions['refitted_total_quadratic']=dict(value=None if beta is None else float(features(s['index'])@beta),
                reason='globally_rank_deficient_no_2D_fit' if beta is None else None)
            for h,r in zip((.01,.02),rawrows):
                v=r['points'][s['id']]
                predictions[f'original_H_{h:.2f}']=dict(value=v['model_energy_change_eV_cell'] if v['passes_original_geometry_guard'] else None,
                    reason=None if v['passes_original_geometry_guard'] else 'model_geometry_outside_original_guard')
            y=observed.get(s['id'])
            for name,p in predictions.items():
                p['raw_error_meV_cell']=None if p['value'] is None or y is None else 1000*(p['value']-y)
                p['shape_error_meV_cell']=None if p['value'] is None or y is None else 1000*((p['value']-origins[name])-(y-observed[0]))
            points.append(dict(id=s['id'],index=s['index'],q_A=s['q_A'],observed=y,models=predictions))
        labelled={p['id'] for p in points if p['observed'] is not None}
        common={p['id'] for p in points if p['id'] in labelled and all(v['value'] is not None for v in p['models'].values())}
        elementary={p['id'] for p in points if p['id'] in labelled and all(p['models'][m]['value'] is not None for m in MODELS if m!='refitted_total_quadratic')}
        proofname,proofpin=PROOFS[channel];proof=pinned(SOURCE.parents[1]/proofname,proofpin)
        costs={r['id']:r['completed_SCFs'] for r in proof['candidate_rows']};costs[0]=0
        result.append(dict(channel=channel,coarse_ids=sorted(coarse),available_coarse_ids=sorted(coarse&set(observed)),
            missing_coarse_ids=sorted(coarse-set(observed)),observed_fine_ids=sorted(labelled),
            original_fine_holes=[p['id'] for p in points if p['observed'] is None],quadratic_design_rank=rank,
            fitted_quadratic_parameter_count=6,quadratic_coefficients=None if beta is None else beta.tolist(),
            model_origin_bias_meV_cell={m:None if v is None else 1000*(v-observed[0]) for m,v in origins.items()},
            common_all_models_ids=sorted(common),common_without_rank_deficient_fit_ids=sorted(elementary),
            secondary_cohort_is_axis_only=all(0 in specs[i]['index'] for i in elementary),points=points,
            common_all_model_statistics={m:statistics(points,common,m) for m in MODELS},
            elementary_common_statistics={m:statistics(points,elementary,m) for m in MODELS if m!='refitted_total_quadratic'},
            original_completed_SCFs_at_selected_coarse_nodes=sum(costs[i] for i in coarse),
            original_SCFs_at_coarse_nodes_with_labels=sum(costs[i] for i in coarse if i in observed),
            original_SCFs_at_missing_coarse_nodes=sum(costs[i] for i in coarse if i not in observed),
            original_completed_SCFs_at_remaining_fine_nodes=sum(costs[p['id']] for p in points),
            original_SCFs_at_fine_holes=sum(costs[p['id']] for p in points if p['observed'] is None),
            source_candidate_SCFs=sum(costs.values()),
            coarse_start_parents=[dict(id=r['id'],parent=r['continuation_parent_id']) for r in proof['candidate_rows'] if r['id'] in coarse],
            coarse_nodes_started_from_fine_parents=[r['id'] for r in proof['candidate_rows']
                if r['id'] in coarse and r['continuation_parent_id'] not in coarse],
            coarse_geometry_and_branch_history_were_not_blinded_to_fine_preparation=True,
            no_2D_surface_claim_from_axis_only_subset=True))
    return dict(status='post_hoc_training_coarsening_not_independent_validation',source_sha256=PIN,
        rows=result,energy_unit='eV per12atomHf4O8cell; error summaries meV/cell',formula_units=4,
        coarse_selection='both original integer grid indices even; no energy/error selection',
        all_labels_already_seen_before_this_exploratory_split=True,models_chosen_after_seen_training=True,
        no_independent_prediction_or_branch_stability_or_physical_error_bound=True,
        interpolation_weights_do_not_certify_geometry_branch_or_error=True,
        fitted_interpolants_and_unfitted_H_models_have_different_information_and_preparation=True,
        no_equal_budget_or_DFT_call_saving_or_B5_novelty_claim=True,
        original296curvature_and118_172face_calls_and_Gamma_canary_G2_preparation_additional=True,
        new_DFT_scheduler_writes_reservations_or_independent_slots_used=0,
        physical_reserved_geometry_generated_or_read=False)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    if args.output.exists():raise FileExistsError('fresh exploratory report required')
    result=analyze()
    with args.output.open('x',encoding='utf-8',newline='\n') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps([dict(channel=r['channel'],coarse_labels=len(r['available_coarse_ids']),
        common_points=len(r['common_all_models_ids']),secondary_points=len(r['common_without_rank_deficient_fit_ids']),
        max_shape_error_meV_cell={m:v['max_absolute_seen_origin_shape_error_meV_cell'] for m,v in r['common_all_model_statistics'].items()}) for r in result['rows']]))
