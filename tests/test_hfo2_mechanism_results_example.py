"""The final mechanism example is offline and keeps local/path targets apart."""
import csv
import json
import pytest
from examples.hfo2_mechanism_results import export, report


def test_missing_evidence_does_not_create_a_plausible_result(tmp_path):
    with pytest.raises(FileNotFoundError):
        export(tmp_path / 'output', root=tmp_path / 'missing_checkout')
    assert not (tmp_path / 'output').exists()


def test_complete_native_report_and_energy_arithmetic(tmp_path):
    data = report()
    assert len(data['local_cross']) == 4
    assert len(data['whole_paths']) == 10
    assert data['replay_new_DFT_calls'] == 0
    assert data['local_points_are_not_new_barriers']
    assert not data['stationary_TS_or_global_ordering_certified']
    assert not data['native_logs_reparsed_by_this_export']
    for row in data['local_cross']:
        assert row['local_candidate_minus_free_sampled_peak_meV_fu'] == pytest.approx(
            row['relative_energy_meV_fu'] - row['connected_free_sampled_peak_meV_fu'])
    out = tmp_path / 'replay'
    assert export(out) == data
    assert json.loads((out / 'mechanism.json').read_bytes()) == data
    with (out / 'local_cross.csv').open(encoding='utf-8', newline='') as handle:
        assert len(list(csv.DictReader(handle))) == 4
    with pytest.raises(FileExistsError):
        export(out)
