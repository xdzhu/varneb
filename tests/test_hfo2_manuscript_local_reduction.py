"""Native-derived local tables; preserve the pre-existing path evidence."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import re

import pytest

ROOT = Path(__file__).resolve().parents[1]
CASES = ROOT / 'benchmarks/hfo2_channels/20261008'
PAPER = ROOT / 'paper/VARNEB_JCTC'
EVIDENCE = CASES / 'manuscript_local_reduction_E162_20261015'
CURVATURE = CASES / 'G3_complete_curvature_E159_20261014'
CANARY = CASES / 'G3_conditional_canary_E160_20261014'


def pinned(path, expected):
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != expected:
        raise ValueError(f'source bytes changed: {path.name}')
    return json.loads(raw)


def sources():
    curvature = pinned(CURVATURE / 'complete_summary.json',
        'f8ca09ed78fbf2102223b66455bbea4e92d84583354d7f792c42a1112e87d25b')
    references = pinned(CURVATURE / 'reference_complements_complete.json',
        '3c5ab39a021ed849ecfa7c1641227ade3f5f21367b272fbb0f135e4562f91f25')
    audits = [pinned(CANARY / f'{channel}_native_audit.json', sha) for channel, sha in (
        ('T', 'b4097af45e0de951c32e550141d745a2d73b6e13eeb4aa8fb3a62be302afa685'),
        ('preserving', '336071744885656b25515f28ca9a6d6b4ecd09f9b411a739827ec7d1e56e5403'))]
    return curvature, references, audits


def expected_rows():
    curvature, references, audits = sources()
    labels = {'T': 'T', 'preserving': 'Preserving'}
    result = []
    for center in curvature['rows']:
        for block in center['blocks']:
            first = ' | '.join(f'{v:.6f}' for v in block['symmetrized_internal_eigenvalues_eV_A2'][:3])
            result.append(f"| {labels[center['channel']]} | {block['step_A']:.2f} | {first} | "
                          f"{block['raw_negative_count']} | {block['tangent_complement_raw_negative_count']} |")
    names = [('T_geometric_triplet', 'T geometric triplet'),
             ('T_two_complete_optical_doublets', 'T complete optical doublets'),
             ('T_all_optical', 'All T optical'),
             *[(f'Cmma_frame_{n}', f'Cmma {n}') for n in (103, 106, 114, 127)]]
    for key, label in names:
        rows = [next(r for r in references['rows'] if
                     (r['reference'], r['channel'], r['step_A']) == (key, center, step))
                for center, step in [('T', .01), ('T', .02), ('preserving', .01), ('preserving', .02)]]
        rank = rows[0]['retained_rank']
        assert all(r['retained_rank'] == rank and r['joint_release']['released_dimension'] == 36-rank for r in rows)
        cells = ' | '.join(f"{r['joint_release']['raw_minimum_eigenvalue_eV_A2']:.6f} "
                           f"({r['joint_release']['raw_negative_count']})" for r in rows)
        result.append(f'| {label} | {rank} | {36-rank} | {cells} |')
    for audit in audits:
        for arm in audit['arms']:
            candidate = arm['final_candidate']
            label = 'Identity, budget exhausted' if candidate is None else 'Curvature seed'
            norm = '-' if candidate is None else f"{candidate['released_gradient_norm_eV_A']:.5f}"
            result.append(f"| {labels[audit['channel']]} | {label} | {arm['completed_SCFs']} | "
                          f"{norm} | {arm['sum_audited_call_seconds']:.2f} |")
    return result


def require_rows(text):
    for row in expected_rows():
        if text.splitlines().count(row) != 1:
            raise ValueError(f'missing, duplicate or changed native-derived row: {row}')


def test_all_fifteen_new_rows_are_native_derived_and_unique():
    assert len(expected_rows()) == 15
    require_rows((PAPER / 'MANUSCRIPT_DRAFT.md').read_text(encoding='utf-8'))


@pytest.mark.parametrize('table', ['| T | 0.01 |', '| All T optical |', '| Cmma 114 |', '| T | Curvature seed |'])
def test_changed_curvature_reference_or_cost_row_fails_closed(table):
    text = (PAPER / 'MANUSCRIPT_DRAFT.md').read_text(encoding='utf-8')
    row = next(r for r in expected_rows() if r.startswith(table))
    with pytest.raises(ValueError, match='native-derived row'):
        require_rows(text.replace(row, row.replace('|', '| CHANGED ', 1)))


def test_old_tables_figures_citations_and_wrapper_are_preserved():
    baseline = json.loads((EVIDENCE / 'baseline_preservation.json').read_text(encoding='utf-8'))
    assert baseline['baseline_git_revision'] == '9e553f22ff9d6dbb24f9b5729251daa824e21d83'
    text = (PAPER / 'MANUSCRIPT_DRAFT.md').read_text(encoding='utf-8')
    counts = Counter(hashlib.sha256(line.strip().encode()).hexdigest()
                     for line in text.splitlines() if line.lstrip().startswith('|'))
    assert not Counter(baseline['original_table_line_sha256_counts']) - counts
    figures = re.findall(r'!\[[^\]]*\]\(([^)]+)\)', text)
    assert figures == baseline['original_figure_targets']
    for name, sha in baseline['original_figure_sha256'].items():
        assert hashlib.sha256((PAPER / name).read_bytes()).hexdigest() == sha
    links = re.findall(r'(?<!!)\[[^\]]*\]\(([^)]+)\)', text)
    assert sorted(set(t for t in links if t.startswith(('https://', 'http://')))) == baseline['external_citation_targets']
    assert all((PAPER / t.split('#')[0]).is_file() for t in links if not t.startswith(('https://', 'http://')))
    # Text checkout may be CRLF on Windows; physical INPUT/KPT stay byte-exact below.
    wrapper = (PAPER / 'varneb_JCTC_draft.tex').read_bytes().replace(b'\r\n', b'\n')
    assert hashlib.sha256(wrapper).hexdigest() == baseline['original_wrapper_sha256']


def test_local_energy_changes_and_censored_cost_are_not_path_barriers():
    _, _, audits = sources()
    text = (PAPER / 'MANUSCRIPT_DRAFT.md').read_text(encoding='utf-8')
    assert sum(a['raw_SCFs_reparsed'] for a in audits) == 31
    for audit in audits:
        identity, seed = audit['arms']
        assert identity['final_candidate'] is None and identity['completed_SCFs'] == 12
        delta = seed['final_candidate']['energy_eV_cell'] - audit['center_source']['energy_eV_cell']
        assert f'{abs(delta):.8f} eV per' in text
        assert not audit['conditional_stability_TS_full_face_or_material_speedup_certified']
    for phrase in ('not a whole-NEB acceleration benchmark', 'Post-relaxation released stability remains',
                   'do not replace any sampled G2 barrier',
                   'a censored cost', '296 curvature SCFs'):
        assert phrase in text
    assert 'whole-path speedup' in (PAPER / 'TERMINOLOGY_LEDGER.md').read_text(encoding='utf-8')


@pytest.mark.parametrize('name,sha', [
    ('INPUT', 'dc6684ffa709bf3d953c647272589b05d844c1293477bd144226d495207cb9df'),
    ('KPT', '92b917107d9df11da28a465cc4900f3574956a057ddf9f7ad125e75c956ec508')])
def test_literal_electronic_snapshots_are_byte_exact(name, sha):
    assert hashlib.sha256((EVIDENCE / name).read_bytes()).hexdigest() == sha
