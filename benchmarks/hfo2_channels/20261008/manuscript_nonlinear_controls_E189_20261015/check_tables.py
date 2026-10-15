"""Bind the manuscript's post-hoc controls to existing scalar evidence."""
from pathlib import Path
import argparse,hashlib,json

ROOT=Path(__file__).resolve().parents[4]
SOURCE=ROOT/'benchmarks/hfo2_channels/20261008/G3_coarsening_E184_20261015/analysis.json'
PIN='377f2c3e2b6c21f0f5d033b2e180ed233fea7547924f2ac5ab993aa4b041f8a5'
LABELS=[('bilinear','Bilinear'),('tensor_quadratic','Tensor quadratic'),
        ('refitted_total_quadratic','Refitted total quadratic'),
        ('original_H_0.01','Original H/g, 0.01 Angstrom'),
        ('original_H_0.02','Original H/g, 0.02 Angstrom')]


def source():
    raw=SOURCE.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=PIN:raise ValueError('original E184 bytes differ')
    return json.loads(raw)


def expected_rows(a):
    t,p=a['rows']
    if (not a['all_labels_already_seen_before_this_exploratory_split']
            or not a['models_chosen_after_seen_training']
            or a['new_DFT_scheduler_writes_reservations_or_independent_slots_used']!=0
            or len(t['available_coarse_ids'])!=9 or len(t['common_all_models_ids'])!=16
            or p['quadratic_design_rank']!=5 or p['quadratic_coefficients'] is not None
            or p['common_all_models_ids']!=[] or not p['secondary_cohort_is_axis_only']
            or any(len(r['coarse_nodes_started_from_fine_parents'])!=8 for r in (t,p))):
        raise ValueError('source does not support the stated retrospective coverage')
    rows=[]
    for name,label in LABELS:
        stat=t['common_all_model_statistics'][name]
        if stat['ids']!=t['common_all_models_ids'] or stat['count']!=16:
            raise ValueError('same sixteen-target cohort required')
        rows.append(f"| {label} | {stat['max_absolute_seen_origin_shape_error_meV_cell']:.6f} | "
                    f"{stat['seen_origin_shape_RMSE_meV_cell']:.6f} |")
    return rows


def require_rows(text,a):
    rows=expected_rows(a)
    if any(text.splitlines().count(r)!=1 for r in rows):
        raise ValueError('source-derived nonlinear row missing, changed or duplicated')
    return rows


def check():
    paper=ROOT/'paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.md'
    a=source();rows=require_rows(paper.read_text(encoding='utf-8'),a)
    return dict(status='posthoc_nonlinear_manuscript_rows_source_checked_not_independent_accuracy',
        source_sha256=PIN,table_rows=rows,new_quantitative_rows=5,
        current_manuscript_sha256=hashlib.sha256(paper.read_bytes()).hexdigest(),
        new_DFT_or_scheduler_writes=0,independent_prediction_gain_or_JCTC_readiness_established=False)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    result=check()
    with a.output.open('x',encoding='utf-8',newline='\n') as s:json.dump(result,s,indent=2);s.write('\n')
    print(json.dumps(result))


if __name__=='__main__':main()
