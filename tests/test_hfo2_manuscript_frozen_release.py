"""Completed matched-ablation manuscript values and explicit evidential limits."""
from pathlib import Path
import json

import pytest

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / 'paper/VARNEB_JCTC'
CASE = ROOT / 'benchmarks/hfo2_channels/20261008/G3_frozen_release_E173_20261015'


def matched_rows(analysis):
    """Complete old-treatment row identities, even when other tables share scalars."""
    expected_rows = []
    for row in analysis['rows']:
        label = 'Preserving' if row['channel'] == 'preserving' else 'T'
        expected = (
            f"| {label} | {row['curvature_probe_step_A']:.2f} | "
            f"{row['max_abs_frozen_shape_error_meV_cell']:.6f} | "
            f"{row['max_abs_released_shape_error_meV_cell']:.6f} | "
            f"{row['max_abs_DFT_release_shape_change_meV_cell']:.6f} |"
        )
        expected_rows.append(expected)
    return expected_rows


def require_matched_rows(text, analysis):
    """Compare entire signed-source-derived table rows, not isolated substrings."""
    for expected in matched_rows(analysis):
        if text.splitlines().count(expected) != 1:
            raise ValueError('one exact source-derived row per center and probe step required')


def test_four_rows_are_actual_both_centers_both_steps_not_selected_maximum_difference():
    analysis = json.loads((CASE / 'analysis.json').read_bytes())
    main = (PAPER / 'MANUSCRIPT_DRAFT.md').read_text(encoding='utf-8')
    require_matched_rows(main, analysis)
    prose = ' '.join(main.split())
    for phrase in (
        'All 16 model-point records are retained',
        'maxima may occur at different coordinates',
        'training-only diagnostics, not independent predictions',
        'not isolate anharmonicity',
        'not the full T/Cmma',
    ):
        # The last two restrictions are expressed without forcing one prose synonym.
        if phrase == 'not isolate anharmonicity':
            assert 'does not isolate anharmonicity' in prose
        elif phrase == 'not the full T/Cmma':
            assert 'Nor does it replace the full T/Cmma' in prose
        else:
            assert phrase in prose
    assert 'twelve-atom cell' in prose and analysis['formula_units'] == 4


@pytest.mark.parametrize('value', ['0.075880','4.671404','4.606959','0.116535',
                                  '4.720758','0.416264','1.327177','1.168761',
                                  '0.403160','1.346071'])
def test_changed_or_duplicated_table_rows_are_rejected(value):
    analysis = json.loads((CASE / 'analysis.json').read_bytes())
    text = (PAPER / 'MANUSCRIPT_DRAFT.md').read_text(encoding='utf-8')
    with pytest.raises(ValueError, match='exact source-derived'):
        require_matched_rows(text.replace(value, '999.999999'), analysis)
    # The new atomic table legitimately reuses F/J scalars. Duplicate the
    # source-bound matched F/J row, not whichever table now occurs first.
    row = next(line for line in matched_rows(analysis) if value in line)
    with pytest.raises(ValueError, match='exact source-derived'):
        require_matched_rows(text + '\n' + row, analysis)


def test_signed_examples_keep_cancellation_and_actual_costs_not_accuracy_claims():
    a = json.loads((CASE / 'analysis.json').read_bytes())
    text = ' '.join((PAPER / 'MANUSCRIPT_DRAFT.md').read_text(encoding='utf-8').split())
    p = a['rows'][0]['points'][0]
    t = a['rows'][2]['points'][3]
    for row in (p, t):
        for key in ('frozen_shape_error_meV_cell', 'missed_release_shape_change_meV_cell',
                    'released_shape_error_meV_cell'):
            assert f'{row[key]:.6f}' in text
    assert 'partly cancelled' in text and 'did not establish general model accuracy' in text
    assert '8.137 completed-call proxy core-hours' in text
    assert '9.076 allocation core-hours' in text
    assert 'curvature and prior candidate preparation remain additional' in text


def test_completed_eight_statics_supersede_preparation_not_the_unfinished_prediction_gate():
    main = ' '.join((PAPER / 'MANUSCRIPT_DRAFT.md').read_text(encoding='utf-8').split())
    methods = ' '.join((PAPER / 'METHODS_DRAFT.md').read_text(encoding='utf-8').split())
    ledger = (PAPER / 'TERMINOLOGY_LEDGER.md').read_text(encoding='utf-8')
    assert 'now completed' in main and 'not completed DFT results' not in main
    assert '22.' in methods and 'eight' in methods and 'supersed' in methods
    assert '| Missed release-shape change |' in ledger
    assert 'independent mechanical-response predictions are still missing' in main
