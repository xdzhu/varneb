"""Check published work diagnostics against frozen material reports, not forecasts."""
import json
from pathlib import Path
import re

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT/'paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.md'
SOURCE = ROOT/'benchmarks/hfo2_channels/20261008/clamped_path_work_E092_20261011'
NAMES = ['zero_M','zero_T','zero_reversing','plus1_T']
LABELS = ['PO-to-M, zero strain','PO-to-T, zero strain',
          'Reversing flip, zero strain','PO-to-T, +1% strain','Preserving flip, zero strain',
          'PO-to-M, +1% strain', 'Preserving flip, +1% strain', 'Reversing flip, +1% strain']
WORK_HEADING = '### E.1. Endpoint closure can conceal a bottleneck-region discrepancy'
NEXT_HEADING = '### E.2. Cross-strain work separates imposed control from configuration change'


def reports():
    old = [json.loads((SOURCE/(name+'.json')).read_text()) for name in NAMES]
    return old + [json.loads((ROOT/'benchmarks/hfo2_channels/20261008'/path).read_text()) for path in
        ('clamped_preserving_reference_E101_20261011/path_work.json',
         'clamped_M_strain_reference_E105_20261012/path_work.json',
         'clamped_plus1_preserving_terminal_E145_20261012/analysis/path_work.json',
         'clamped_plus1_reversing_terminal_E148_20261014/analysis/path_work.json')]


def section(text=None):
    if text is None:
        text = PAPER.read_text(encoding='utf-8')
    assert text.count(WORK_HEADING) == text.count(NEXT_HEADING) == 1
    body = text.split(WORK_HEADING, 1)[1].split(NEXT_HEADING, 1)[0]
    assert not re.search(r'^#{1,3} ', body, re.M)
    return body


def cmma_ranges(text=None):
    passage = section(text).split('ranges are ', 1)[1].split('Each representation', 1)[0]
    return [float(value) for value in re.findall(r'\d+\.\d{6}', passage)]


@pytest.mark.parametrize('index',range(8))
def test_table7_peak_and_endpoint_values_keep_sign_order_and_formula_unit(index):
    data = reports()[index]
    text = section()
    row = next(line for line in text.splitlines() if line.startswith('| '+LABELS[index]+' |'))
    values = [float(value.strip()) for value in row.split('|')[2:-1]]
    peak = data['sampled_peak_image']
    span = data['representations']['T_geometric_triplet']
    factor = 1000/data['formula_units']
    assert factor == 250
    expected = factor*np.asarray([span['cumulative_energy_change'][peak],
        span['cumulative_full_work'][peak],span['cumulative_energy_minus_work'][peak],
        span['endpoint_energy_minus_work']])
    np.testing.assert_allclose(values,expected,atol=5.0001e-7,rtol=0)
    assert data['new_DFT_calls'] == data['scheduler_mutations'] == 0
    assert data['G3_selection_or_holdout_access'] is False


@pytest.mark.parametrize('marker,end,space,part', [
    ('signed retained work to the respective peaks is','meV/HfO2','T_geometric_triplet','retained_atomic'),
    ('Omitted atomic work is','; open-cell','T_geometric_triplet','omitted_atomic'),
    ('open-cell work is','The explicit','T_geometric_triplet','open_cell'),
    ('doublets yield retained work','in the same order','T_two_lowest_complete_optical_doublets','retained_atomic'),
])
def test_signed_work_prose_matches_complete_declared_span(marker,end,space,part):
    passage = section().split(marker,1)[1].split(end,1)[0]
    actual = [float(value) for value in re.findall(r'-?\d+\.\d{6}',passage)]
    expected = [250*data['representations'][space]['cumulative_partition_work'][part][data['sampled_peak_image']]
                for data in reports()]
    np.testing.assert_allclose(actual,expected,atol=5.0001e-7,rtol=0)
    if space == 'T_geometric_triplet' and part == 'retained_atomic':
        assert actual[0] < 0  # signed work is not truncated to an energy fraction


def test_all_four_Cmma_registrations_are_kept_as_ranges_not_selected_winners():
    actual = cmma_ranges()
    expected = []
    for data in reports():
        span = data['representations']
        values = [250*span['Cmma_frame_'+frame]['cumulative_partition_work']['retained_atomic'][data['sampled_peak_image']]
                  for frame in ('103','106','114','127')]
        expected.extend((min(values),max(values)))
        assert span['T_geometric_triplet']['retained_rank'] == 3
        assert all(span['Cmma_frame_'+frame]['retained_rank'] == 4 for frame in ('103','106','114','127'))
    np.testing.assert_allclose(actual,expected,atol=5.0001e-7,rtol=0)


def test_native_chart_agreement_is_only_an_algebraic_check_and_reused_frames_not_replicates():
    rows = reports()
    assert sum(row['portable_native_frames_reparsed'] for row in rows) == 72
    maximum = max(span['independent_native_work_max_defect_eV_cell']
                  for row in rows for span in row['representations'].values())
    assert maximum <= 1.12e-16
    translation_peak_max = max(abs(250*span['cumulative_partition_work']['rigid_translation'][row['sampled_peak_image']])
                  for row in rows for span in row['representations'].values())
    assert translation_peak_max <= 4.63e-9
    assert all(span['physical_error_bound_measured'] is False
               for row in rows for span in row['representations'].values())
    assert '66 distinct raw-log/geometry' in section()
    assert 'algebraic agreement does not test DFT accuracy' in section()


def test_section_numbers_and_new_evidence_links_are_unambiguous():
    text = PAPER.read_text(encoding='utf-8')
    numbers = re.findall(r'^### ((?:\d+|[A-Z])\.\d+)\.',text,re.M)
    assert len(numbers) == len(set(numbers))
    # The mechanism-first manuscript moved the former work-theory methods
    # into G.3/G.4; retain their content/link checks, not obsolete numbering.
    assert [n for n in numbers if n.startswith('2.')] == ['2.1', '2.2', '2.3', '2.4']
    assert all(n in numbers for n in ('G.3', 'G.4', 'E.1', 'E.2'))
    assert [n for n in numbers if n.startswith('3.')] == ['3.1', '3.2', '3.3', '3.4']
    assert 'Table A11' in section()
    assert '### G.3. Signed path work as an energy-consistency diagnostic' in text
    assert '### G.4. Work between the two seen strain conditions' in text
    method = text.split('### G.3. ',1)[1].split('### G.5.',1)[0]
    for target in re.findall(r'\]\(([^)]+)\)',method+section()):
        assert not target.startswith(('http','file:'))
        assert (PAPER.parent/target).resolve().is_file()


@pytest.mark.parametrize('heading', [WORK_HEADING, NEXT_HEADING])
@pytest.mark.parametrize('fault', ['missing', 'duplicate'])
def test_work_section_requires_unique_current_boundaries(heading, fault):
    text = PAPER.read_text(encoding='utf-8')
    broken = text.replace(heading, '') if fault == 'missing' else text+'\n'+heading
    with pytest.raises(AssertionError):
        section(broken)


def test_work_table_cannot_borrow_a_missing_row_from_another_section():
    text = PAPER.read_text(encoding='utf-8')
    body = section(text)
    row = next(line for line in body.splitlines() if line.startswith('| '+LABELS[0]+' |'))
    broken = text.replace(body, body.replace(row, ''), 1)+'\n'+row
    with pytest.raises(StopIteration):
        next(line for line in section(broken).splitlines()
             if line.startswith('| '+LABELS[0]+' |'))


def test_each_Cmma_bound_still_rejects_mutation_despite_correct_values_elsewhere():
    text = PAPER.read_text(encoding='utf-8')
    body = section(text)
    expected = cmma_ranges(text)
    assert len(expected) == 16
    for value in expected:
        changed = body.replace(f'{value:.6f}', f'{value+1:.6f}', 1)
        broken = text.replace(body, changed, 1)+'\n'+body
        with pytest.raises(AssertionError):
            np.testing.assert_allclose(cmma_ranges(broken), expected, atol=5.0001e-7, rtol=0)
