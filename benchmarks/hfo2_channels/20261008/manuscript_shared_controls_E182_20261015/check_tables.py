"""Bind manuscript shared-observable tables to the untouched E180 scalars."""
from pathlib import Path
import argparse,hashlib,json,sys
ROOT=Path(__file__).resolve().parents[4]
SOURCE=ROOT/'benchmarks/hfo2_channels/20261008/G3_shared_Q_controls_E180_20261015/analysis.json'
SOURCE_SHA='39e31d790592e7a2f640a8428f00efecc28fd14775b7fb4b2b699b7da54c6402'
LABELS={'T_geometric_triplet':'T geometric triplet',
    'T_two_complete_optical_doublets':'T complete optical doublets','T_all_optical':'All T optical',
    **{f'Cmma_frame_{i}':f'Cmma {i}' for i in (103,106,114,127)},
    'adaptive_two_directions_plus_complete_complement':'Coupled full-space'}


def source():
    raw=SOURCE.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=SOURCE_SHA:raise ValueError('original E180 scalar bytes differ')
    return json.loads(raw)


def expected_rows(a):
    if (a['model_point_rows']!=800 or a['fitted_coefficients']!=0 or a['new_DFT_or_scheduler_writes']!=0
            or a['stable_conditional_branch_TS_physical_error_or_independent_advantage_certified']):
        raise ValueError('source does not support the reported evidence scope')
    rows=[]
    for r in a['rows']:
        if r['reference']=='T_all_optical':
            stats=r['own_geometry_eligible_statistics'];name='T' if r['channel']=='T' else 'Preserving'
            rows.append(f"| {name} | {r['step_A']:.2f} | {stats['count']}/{stats['denominator_non_origin']} | "
                f"{stats['max_absolute_seen_origin_aligned_shape_discrepancy_meV_cell']:.6f} |")
    for reference,label in LABELS.items():
        selected=[r for r in a['rows'] if r['reference']==reference]
        if len(selected)!=4:raise ValueError('both centers and both steps required')
        if any(r['common_all_reference_geometry_eligible_statistics']['count']!=0 for r in selected):
            raise ValueError('claimed empty common nonzero set differs')
        if any(r['common_all_reference_geometry_eligible_statistics']['max_absolute_energy_discrepancy_meV_cell'] is not None for r in selected):
            raise ValueError('missing statistics must remain null, not zero')
        percenter={}
        for channel in ('T','preserving'):
            pair=[r for r in selected if r['channel']==channel]
            if len(pair)!=2 or {r['step_A'] for r in pair}!={.01,.02}:raise ValueError('two steps required')
            values=[r['own_geometry_eligible_statistics']['count'] for r in pair]
            if values[0]!=values[1]:raise ValueError('step-independent eligibility count no longer holds')
            percenter[channel]=values[0]
        r=selected[0];rank='-' if r['nominal_atomic_reference_rank'] is None else str(r['nominal_atomic_reference_rank'])
        rows.append(f"| {label} | {rank} | {r['allowed_dimension']} | {percenter['T']}/24 | {percenter['preserving']}/24 |")
    return rows


def require_rows(text,a):
    rows=expected_rows(a)
    if len(rows)!=12:raise ValueError('exact twelve new quantitative rows required')
    for row in rows:
        if text.splitlines().count(row)!=1:raise ValueError('source-derived row missing, changed or duplicated')
    return rows


def check():
    a=source();paper=ROOT/'paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.md'
    rows=require_rows(paper.read_text(encoding='utf-8'),a)
    return dict(status='manuscript_tables_bound_to_seen_shared_observable_controls',
        source_sha256=SOURCE_SHA,new_table_rows=rows,new_row_count=len(rows),
        source_model_records=800,new_DFT_or_scheduler_writes=0,
        complete_reference_energy_equivalence_eV_cell=max(r['max_complete_T_vs_adaptive_Schur_energy_difference_eV_cell'] for r in a['complete_space_equivalence']),
        complete_reference_coordinate_equivalence_A=max(r['max_complete_T_vs_adaptive_Schur_coordinate_difference_A'] for r in a['complete_space_equivalence']),
        independent_accuracy_stability_or_JCTC_readiness_established=False,
        current_manuscript_sha256=hashlib.sha256(paper.read_bytes()).hexdigest())


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    result=check()
    with a.output.open('x',encoding='utf-8',newline='\n') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result))


if __name__=='__main__':main()
