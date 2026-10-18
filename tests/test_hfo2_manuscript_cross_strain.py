"""Paper tables remain native-derived after narrative reordering, not forecasts."""
import json
from pathlib import Path
import re

import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
PAPER=ROOT/'paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.md'
SOURCE=ROOT/'benchmarks/hfo2_channels/20261008/cross_strain_work_E151_20261014/analysis.json'
LABELS={'M':'PO-to-M','T':'PO-to-T','preserving':'Preserving flip','reversing':'Reversing flip'}


@pytest.mark.parametrize('candidate',list(LABELS))
def test_cross_strain_table_has_actual_signed_work_and_unit_normalization(candidate):
    data=json.loads(SOURCE.read_text())
    row=next(r for r in data['candidates'] if r['candidate']==candidate)
    pair=row['matched_peak_pair']; parts=pair['relative_PO_parts_meV_fu']
    text=PAPER.read_text(encoding='utf-8').split('### E.2.',1)[1].split('### E.3.',1)[0]
    line=next(line for line in text.splitlines() if line.startswith('| '+LABELS[candidate]+' |'))
    actual=[float(s.strip()) for s in line.split('|')[2:-1]]
    expected=[row['observed_sampled_barrier_change_meV_fu'],parts['internal_atomic'],
        parts['open_cell'],parts['prescribed_substrate'],pair['relative_PO_energy_minus_work_meV_fu']]
    np.testing.assert_allclose(actual,expected,atol=5.0001e-7,rtol=0)
    assert pair['relative_PO_energy_change_meV_fu']==pytest.approx(expected[0],abs=1e-12)
    assert abs(parts['rigid_translation'])<2.411e-10


def test_results_lead_with_whole_matrix_and_keep_all_original_nine_figures():
    text=PAPER.read_text(encoding='utf-8')
    main=text.split('## 3.',1)[1].split('## 4.',1)[0]
    assert re.findall(r'^### (3\.\d+)\.',main,re.M)==['3.1','3.2','3.3','3.4']
    assert main.index('Strain changes channel competition')<main.index('Distinct oxygen-shuffle routes')
    assert '0.464068' in main and 'not the favorable signature' in main
    assert 'not an error-certified rejection' in main
    assert len(re.findall(r'^!\[',text,re.M))==14
    labels=re.findall(r'^Figure (\w+)\.',text,re.M)
    assert labels==['1','2','3','4','5','6','A1','A2','A3','A4','A5','A6','A7','A8']
    assert 'figures/hfo2_matched_T_20261018/matched_T_release.pdf' in main
    for target in re.findall(r'^!\[[^\]]*\]\(([^)]+)\)',text,re.M):
        assert (PAPER.parent/target).is_file()
    # Detailed barrier/coverage tables are relocated, never dropped.
    for name in ('Table A1.','Table A2.','Table A3.','Table A4.','Table A5.'):
        assert name in text
    assert '## Appendix C.' in text and '## Appendix D.' in text


def test_cross_strain_work_and_missing_prediction_are_not_promoted_to_mechanisms():
    text=PAPER.read_text(encoding='utf-8')
    cross=text.split('### E.2.',1)[1].split('### E.3.',1)[0]
    for phrase in ('not a barrier-error estimate','not a\nmeasurement that a frozen strain slice',
        'No intermediate\ncross-strain geometries were calculated','not a total DFT'):
        assert phrase in (text if phrase=='not a total DFT' else cross)
    discussion=text.split('## 4.',1)[1].split('## Appendix A.',1)[0]
    normalized=' '.join(discussion.split())
    assert 'A complete change of coordinates leaves the same quadratic response unchanged' in normalized
    assert 'A new predictive-method claim would require additional independent evidence' in normalized
