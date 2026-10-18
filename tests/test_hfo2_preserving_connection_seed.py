"""The connection seed must not pin its center or create new DFT work."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest
from ase.io import read
from examples.hfo2_fixed_input_factory import CONTRACT, FixedHfo2Calculator
from scripts.audit_hfo2_static_replica import sha256

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT/'benchmarks/hfo2_channels/20261008/preserving_connection_E273_20261019'
spec = importlib.util.spec_from_file_location('connection_seed_test', CASE/'prepare.py')
preparer = importlib.util.module_from_spec(spec); spec.loader.exec_module(preparer)


@pytest.fixture(scope='module')
def prepared(tmp_path_factory):
    output = tmp_path_factory.mktemp('P_connection')/'prepared'
    return output, preparer.prepare(output)


def test_one_chain_has_two_fixed_endpoints_and_seven_free_interiors(prepared):
    output, meta = prepared; images = read(output/'seed.traj', ':')
    assert len(images) == 9 and meta['n_fixed_endpoints'] == 2
    assert meta['free_interiors'] == list(range(1, 8))
    assert meta['retained_pattern_constraints_in_NEB'] == []
    assert meta['fmax_target_eV_A'] == .10 and not meta['climb']
    assert meta['maximum_updates'] == 24 and meta['maximum_new_interior_SCFs'] == 168
    assert meta['walltime'] == '06:00:00' and meta['maximum_study_allocations'] == 2
    for name, image in (('initial.vasp', images[0]), ('final.vasp', images[-1])):
        other = read(output/name)
        np.testing.assert_allclose(other.positions, image.positions, atol=2e-13, rtol=0)
        np.testing.assert_allclose(other.cell.array, image.cell.array, atol=2e-13, rtol=0)


def test_nine_exact_caches_and_same_six_physical_inputs(prepared):
    output, _ = prepared
    parameters = json.loads((output/'factory_parameters.json').read_bytes())
    assert set(parameters) == {'source_directory', 'seed_cache_records'}
    records = parameters['seed_cache_records']; assert len(records) == 9 and all(records)
    for i, record in enumerate(records):
        assert set(record) == {'directory', 'input_sha256', 'raw_log_sha256'}
        assert all(record['input_sha256'][n] == h for n, h in CONTRACT.items())
        if i == 4:
            assert '/results/0/minus/evaluations/eval-000014-' in record['directory']
        else:
            assert 'preserving-branches-E272' not in record['directory']


def test_preparation_is_not_execution_and_all_generated_bytes_match(prepared):
    output, meta = prepared
    assert meta['new_DFT_calls'] == meta['scheduler_mutations'] == 0
    assert not meta['execution_admitted'] and meta['source_job_must_finish'] == '28948425'
    assert meta['original_labels_and_scores_unchanged']
    for name, digest in meta['files_sha256'].items():
        assert sha256(output/name) == digest
    assert meta['preparation_source_sha256'] == sha256(CASE/'prepare.py')
    assert meta['protocol_sha256'] == sha256(CASE/'protocol.md')
    with pytest.raises(FileExistsError):
        preparer.prepare(output)


def test_preparation_cannot_invoke_native_DFT(tmp_path, monkeypatch):
    monkeypatch.setattr(FixedHfo2Calculator, 'calculate', lambda *a, **kw: pytest.fail('no DFT in preparation'))
    meta = preparer.prepare(tmp_path/'new')
    assert meta['new_DFT_calls'] == 0 and not meta['execution_admitted']
