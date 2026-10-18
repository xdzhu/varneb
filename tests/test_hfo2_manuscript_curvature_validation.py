"""Cross-check manuscript tables against the unchanged E164/E165 analyses."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import re

import pytest

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'benchmarks/hfo2_channels/20261008'
PAPER=ROOT/'paper/VARNEB_JCTC'
CASE=BASE/'manuscript_curvature_validation_E167_20261015'


def pinned(relative,expected):
    raw=(BASE/relative).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=expected:
        raise ValueError('original analysis bytes changed')
    return json.loads(raw)


def sources():
    return (pinned('G3_atomic_cell_signature_E164_20261015/analysis.json',
                   'ef7798df3e71fbbc3e86bd468b608e4bf732cb13cd048e1518b8608d4753db3a'),
            pinned('G3_quadratic_agreement_E165_20261015/analysis.json',
                   'c2e51c37b185a00676a6ff21660ea132384c0a22cfac2cebd18e30a5e5a9b745'))


def expected_rows():
    atomic,models=sources();rows=[]
    labels={'T':'T','preserving':'Preserving'}
    for r in atomic['rows']:
        rows.append(f"| {labels[r['channel']]} | {r['step_A']:.2f} | "
            f"{r['atomic_fixed_cell_eigenvalues_eV_A2'][0]:.6f} | "
            f"{r['internal_eigenvalues_eV_A2'][0]:.6f} | "
            f"{r['open_cell_eigenvalues_eV_A2'][0]:.6f} | "
            f"{r['negative_subspace']['open_cell_projector_trace']:.6g} |")
    for r in models['rows']:
        rows.append(f"| {labels[r['channel']]} | {r['curvature_probe_step_A']:.2f} | "
            f"{r['observed_DFT_candidate_origin_release_meV_cell']:.6f} | "
            f"{r['raw_formal_quadratic_origin_release_meV_cell']:.6f} | "
            f"{r['common_numeric_Q_statistics']['max_absolute_zero_aligned_shape_error_meV_cell']:.6f} |")
    return rows


def require_rows(text):
    for row in expected_rows():
        if text.splitlines().count(row)!=1:
            raise ValueError('missing, repeated or changed source-derived row')


def test_eight_new_rows_keep_both_original_centers_and_both_steps():
    assert len(expected_rows())==8
    require_rows((PAPER/'MANUSCRIPT_DRAFT.md').read_text(encoding='utf-8'))


@pytest.mark.parametrize('index',range(8))
def test_any_changed_or_repeated_quantitative_row_is_rejected(index):
    text=(PAPER/'MANUSCRIPT_DRAFT.md').read_text(encoding='utf-8')
    row=expected_rows()[index]
    for bad in (text.replace(row,''),text.replace(row,row+'\n'+row),
                text.replace(row,row.replace('|','| altered ',1))):
        with pytest.raises(ValueError,match='source-derived row'):require_rows(bad)


def test_existing_tables_figures_citations_and_wrapper_are_not_displaced():
    baseline=json.loads((CASE/'baseline_preservation.json').read_bytes())
    assert baseline['baseline_git_revision']=='2936fcb03fa4925960755c598016ceec731bba26'
    assert baseline['original_table_line_count']==101
    text=(PAPER/'MANUSCRIPT_DRAFT.md').read_text(encoding='utf-8')
    counts=Counter(hashlib.sha256(line.strip().encode()).hexdigest()
                   for line in text.splitlines() if line.lstrip().startswith('|'))
    assert not Counter(baseline['original_table_line_sha256_counts'])-counts
    figures=re.findall(r'!\[[^\]]*\]\(([^)]+)\)',text)
    # All original assets remain; the matched-T control adds one main figure.
    assert not Counter(baseline['original_figure_targets'])-Counter(figures)
    assert len(figures)==14 and len(set(figures))==14
    assert 'figures/hfo2_matched_T_20261018/matched_T_release.pdf' in figures
    assert all(hashlib.sha256((PAPER/p).read_bytes()).hexdigest()==h
               for p,h in baseline['original_figure_sha256'].items())
    links=re.findall(r'(?<!!)\[[^\]]*\]\(([^)]+)\)',text)
    assert set(baseline['external_citation_targets']) <= {p for p in links if p.startswith('http')}
    assert all((PAPER/p.split('#')[0]).is_file() for p in links if not p.startswith('http'))
    wrapper=(PAPER/'varneb_JCTC_draft.tex').read_bytes().replace(b'\r\n',b'\n')
    # The author-review header can evolve; the single-source rendering stays.
    assert b'\\markdownInput{MANUSCRIPT_DRAFT.md}' in wrapper
    assert b'outputDir={build}' in wrapper and b'\\FloatBarrier' in wrapper


def test_raw_atomic_cell_observations_are_not_a_strain_response_or_TS_claim():
    atomic,_=sources();text=(PAPER/'MANUSCRIPT_DRAFT.md').read_text(encoding='utf-8')
    assert len(atomic['rows'])==4
    for r in atomic['rows']:
        assert r['atomic_raw_negative_count']==r['internal_raw_negative_count']==2
        assert r['open_cell_raw_nonpositive_count']==0
        assert not r['phonon_identity_stationarity_branch_stability_or_barrier_decomposition_certified']
    for phrase in ('not a general rule excluding coupling-',
                   'not arbitrary changes of the cell length scale',
                   'substrate response and the permitted open-cell fluctuations are',
                   'not a percentage of a barrier energy'):
        assert phrase in text


def test_shape_errors_are_signed_snapshot_training_diagnostics_not_independent_accuracy():
    _,models=sources();text=(PAPER/'MANUSCRIPT_DRAFT.md').read_text(encoding='utf-8')
    assert models['common_numeric_Q_ids']==list(range(6))
    assert models['available_nodes']=={'T':9,'preserving':6}
    for r in models['rows']:
        error=next(p['zero_aligned_shape_error_meV_cell'] for p in r['points'] if p['id']==5)
        assert (error>0)==(r['channel']=='T')
        assert r['common_numeric_Q_statistics']['n']==6
        if r['curvature_probe_step_A']==.01:
            assert r['formal_Schur_softening_Frobenius_norm_eV_A2']<1e-24
    for phrase in ('both original probe steps', 'no coefficients',
                   'not a calibration-free or independent prediction',
                   'seen-training diagnostics, not independent errors',
                   'training-only diagnostics, not independent predictions', 'not infinitesimal Hessians',
                   'eight', 'audited source snapshot'):
        assert phrase.lower() in text.lower()
    assert not models['conditional_stability_full_faces_TS_or_material_prediction_gain_certified']
    assert not models['physical_reserved_geometry_generated_or_read']


def test_new_terms_and_methods_record_keep_budget_and_metric_boundaries():
    ledger=(PAPER/'TERMINOLOGY_LEDGER.md').read_text(encoding='utf-8')
    methods=(PAPER/'METHODS_DRAFT.md').read_text(encoding='utf-8')
    for term in ('Open-cell trace','Origin-aligned shape discrepancy','Matched frozen-R control'):
        assert f'| {term} |' in ledger
    for phrase in ('leaving 52 uncommitted of 81','eight independent slots remain untouched',
                   'not a fully','not computed','nonzero cell gradients'):
        assert phrase.lower() in methods.lower()
