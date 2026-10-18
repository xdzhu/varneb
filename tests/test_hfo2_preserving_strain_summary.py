"""Main-figure completeness gate uses the actual partial archive, never mock labels."""
import importlib.util
import json
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT/'benchmarks/hfo2_channels/20261008/preserving_strain_E275_20261019'


def test_partial_actual_result_never_becomes_a_complete_three_strain_figure(tmp_path):
    spec = importlib.util.spec_from_file_location('E275_summary_test', CASE/'plot_summary.py')
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    source = CASE/'completed_analysis_0050/analysis.json'
    data = json.loads(source.read_bytes())
    with pytest.raises(ValueError, match='complete original three-strain'): m.summary_rows(data)
    with pytest.raises(ValueError, match='complete original three-strain'): m.draw(source, tmp_path/'figure')
    assert not (tmp_path/'figure').exists()
