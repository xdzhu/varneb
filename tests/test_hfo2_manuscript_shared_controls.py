"""New main/appendix quantitative rows and explicit no-advantage boundaries."""
from pathlib import Path
from copy import deepcopy
import importlib.util
import pytest

ROOT=Path(__file__).resolve().parents[1]
CASE=ROOT/'benchmarks/hfo2_channels/20261008/manuscript_shared_controls_E182_20261015'
s=importlib.util.spec_from_file_location('E182_table_checks',CASE/'check_tables.py')
a=importlib.util.module_from_spec(s);s.loader.exec_module(a)
PAPER=ROOT/'paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.md'


def test_all_twelve_quantitative_rows_keep_strong_complete_reference_and_all_frames():
    rows=a.require_rows(PAPER.read_text(encoding='utf-8'),a.source())
    assert len(rows)==12 and sum('Cmma' in row for row in rows)==4


@pytest.mark.parametrize('index',range(12))
def test_every_changed_missing_or_duplicated_new_row_is_rejected(index):
    text=PAPER.read_text(encoding='utf-8');data=a.source();row=a.expected_rows(data)[index]
    for broken in (text.replace(row,''),text.replace(row,row+'\n'+row),text.replace(row,row.replace('|','| changed',1))):
        with pytest.raises(ValueError,match='source-derived'):a.require_rows(broken,data)


@pytest.mark.parametrize('fault',['zero_error','common_set','extra_DFT','fitted','false_certification','missing_step'])
def test_overclaims_and_missing_source_dimensions_rejected(fault):
    data=deepcopy(a.source())
    if fault=='zero_error':data['rows'][0]['common_all_reference_geometry_eligible_statistics']['max_absolute_energy_discrepancy_meV_cell']=0.
    elif fault=='common_set':data['rows'][0]['common_all_reference_geometry_eligible_statistics']['count']=1
    elif fault=='extra_DFT':data['new_DFT_or_scheduler_writes']=800
    elif fault=='fitted':data['fitted_coefficients']=1
    elif fault=='false_certification':data['stable_conditional_branch_TS_physical_error_or_independent_advantage_certified']=True
    else:data['rows'].pop()
    with pytest.raises(ValueError):a.expected_rows(data)


def test_prose_retains_equal_problem_null_vs_finite_domain_and_missing_evidence():
    text=' '.join(PAPER.read_text(encoding='utf-8').split())
    for phrase in ('not competing accuracy results','not DFT accuracy',
                   'not replaced with their last optimization trials',
                   'not a new theorem or a source of prediction gain',
                   'not evidence that a native Cmma-coordinate energy model fails',
                   'equal preparation cost or speedup was not established',
                   'independent mechanical-response predictions are still missing',
                   'not a measured physical stability-resolution floor'):
        assert phrase in text
    assert '2.64' in text and '1.56' in text
    ledger=(PAPER.parent/'TERMINOLOGY_LEDGER.md').read_text(encoding='utf-8')
    for term in ('Shared-observable control','Complete-space equivalence','Geometry-eligible comparison set'):
        assert f'| {term} |' in ledger
