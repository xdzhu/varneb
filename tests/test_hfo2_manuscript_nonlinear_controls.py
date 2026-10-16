"""The new comparison must remain retrospective and retain all five arms."""
from pathlib import Path
from copy import deepcopy
from collections import Counter
import hashlib,importlib.util,json,re

import pytest

ROOT=Path(__file__).resolve().parents[1]
CASE=ROOT/'benchmarks/hfo2_channels/20261008/manuscript_nonlinear_controls_E189_20261015'
spec=importlib.util.spec_from_file_location('E189_manuscript_controls',CASE/'check_tables.py')
a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
PAPER=ROOT/'paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.md'


def test_five_rows_use_same_sixteen_points_and_both_unfitted_steps():
    rows=a.require_rows(PAPER.read_text(encoding='utf-8'),a.source())
    assert len(rows)==5
    assert any('Bilinear' in r for r in rows) and sum('Original H/g' in r for r in rows)==2


@pytest.mark.parametrize('index',range(5))
def test_changed_missing_or_duplicate_quantitative_row_is_rejected(index):
    text=PAPER.read_text(encoding='utf-8');data=a.source();row=a.expected_rows(data)[index]
    for broken in (text.replace(row,''),text.replace(row,row+'\n'+row),text.replace(row,row.replace('|','| wrong',1))):
        with pytest.raises(ValueError,match='source-derived'):a.require_rows(broken,data)


@pytest.mark.parametrize('fault',['blind','no_missing','changed_cohort','DFT','no_branch_leakage'])
def test_source_cannot_silently_support_a_stronger_claim(fault):
    data=deepcopy(a.source())
    if fault=='blind':data['all_labels_already_seen_before_this_exploratory_split']=False
    elif fault=='no_missing':data['rows'][1]['quadratic_design_rank']=6
    elif fault=='changed_cohort':data['rows'][0]['common_all_model_statistics']['tensor_quadratic']['count']=15
    elif fault=='DFT':data['new_DFT_scheduler_writes_reservations_or_independent_slots_used']=16
    else:data['rows'][0]['coarse_nodes_started_from_fine_parents']=[]
    with pytest.raises(ValueError):a.expected_rows(data)


def test_prose_retains_bias_coverage_cost_and_missing_independent_evidence():
    text=' '.join(PAPER.read_text(encoding='utf-8').split())
    for term in ('6.120672','+1.546418','not an independent forecast or a saving of DFT calls',
                 'leaving the six-term quadratic design rank five','47/71','58/114',
                 'withholding fine energies from a fit did not blind the preceding workflow',
                 'the smaller quadratic maxima do not represent a model improvement',
                 'not itself a method advance','independent material predictions are still missing'):
        assert term in text
    ledger=(PAPER.parent/'TERMINOLOGY_LEDGER.md').read_text(encoding='utf-8')
    for name in ('Post-hoc coarsening','Tensor-quadratic interpolation'):assert f'| {name} |' in ledger


def test_prior135_table_lines_nine_figures_external_citations_and_wrapper_preserved():
    old=json.loads((CASE/'preservation.json').read_bytes())
    text=PAPER.read_text(encoding='utf-8')
    current=Counter(l for l in text.splitlines() if l.startswith('|'))
    assert len(old['prior_table_lines'])==135
    # User-approved main/appendix reordering changes order, not row values or multiplicity.
    assert not Counter(old['prior_table_lines'])-current
    # All prior rows remain; atomic treatment and four statics add six lines each.
    assert len([l for l in text.splitlines() if l.startswith('|')])==154
    assert 'Table A17. Four frozen directional training statics' in text
    assert re.findall(r'!\[[^\]]*\]\(([^)]+)\)',text)==old['prior_figure_targets']
    assert re.findall(r'\]\((https?://[^)]+)\)',text)==old['prior_external_links']
    for name,sha in old['prior_figure_sha256'].items():
        assert hashlib.sha256((PAPER.parent/name).read_bytes()).hexdigest()==sha
    assert hashlib.sha256((PAPER.parent/'varneb_JCTC_draft.tex').read_bytes()).hexdigest()==old['wrapper_sha256']
