"""The main argument must reflect the two real completed follow-ups."""
import json
from pathlib import Path
import re

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'benchmarks/hfo2_channels/20261008'
PAPER=ROOT/'paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.md'


def test_whole_chain_lowering_and_initialization_difference_are_in_main_results():
    text=PAPER.read_text(encoding='utf-8')
    main=text.split('## 3.',1)[1].split('## 4.',1)[0]
    connection=json.loads((BASE/'preserving_connection_E273_20261019/figures/terminal_step_0009/analysis.json').read_bytes())
    midpoint=json.loads((BASE/'midpoint_reverse_E269_20261018/terminal_comparison/comparison.json').read_bytes())
    assert f"{connection['bands'][-1]['sampled_polygon_peak_meV_fu']:.6f}" in main
    assert f"{midpoint['paired']['per_channel_upper_minus_zero_meV_fu']['preserving']:.6f}" in main
    assert f"{midpoint['paired']['from_upper_T_minus_preserving_meV_fu']:.6f}" in main
    assert 'not resolved as a physical ordering' in main
    assert 'does not survive this initialization control' in main
    assert 'corresponding preserving control is not yet complete' not in text


def test_no_frozen_target_replacement_or_stationary_saddle_claim():
    text=PAPER.read_text(encoding='utf-8')
    assert 'registered training labels are not replaced' in text
    assert 'Both original prospective labels and all frozen forecast scores remain' in text
    assert 'Neither a global MEP nor a stationary' in text
    assert 'the small T releases cannot fill that gap' in text
    figures=re.findall(r'^Figure (\d+)\.',text,re.M)
    assert figures==[str(i) for i in range(1,9)]


def test_frozen_crossover_table_uses_the_measured_conditional_contrasts():
    text=PAPER.read_text(encoding='utf-8')
    data=json.loads((BASE/'midpoint_factorial_E274_20261019/analysis/analysis.json').read_bytes())
    table=text.split('Table 4. Frozen midpoint contrasts',1)[1].split('![',1)[0]
    for label,row in zip(('Preserving','PO-to-T'),data['rows']):
        line=next(l for l in table.splitlines() if l.startswith('| '+label+' |'))
        measured=[float(s.strip()) for s in line.split('|')[2:-1]]
        expected=[row[k] for k in ('diagonal_change_meV_fu','internal_at_A_cell_meV_fu',
            'cell_at_A_internal_meV_fu','finite_interaction_meV_fu')]
        assert len(measured)==len(expected)==4
        assert all(abs(a-b)<=5.0001e-7 for a,b in zip(measured,expected))
    assert 'neither reoptimized barriers nor unique' in table
