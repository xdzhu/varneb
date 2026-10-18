"""Real native source/endpoint geometry and finite runtime, no DFT/scheduler."""
import importlib.util
from pathlib import Path
import json
import numpy as np
import pytest
from ase.io import read
from examples.hfo2_fixed_input_factory import CONTRACT
from vcneb import validate_periodic_path_lift

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT/'benchmarks/hfo2_channels/20261008/preserving_strain_E275_20261019'


def load(name):
    spec = importlib.util.spec_from_file_location('E275_'+name, CASE/(name+'.py'))
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module); return module


@pytest.fixture(scope='module')
def prepared(tmp_path_factory):
    output = tmp_path_factory.mktemp('E275')/'seeds'
    return output, load('prepare').prepare(output)


def test_same_native_source_both_strains_and_original_endpoint_caches(prepared):
    output, report = prepared
    assert report['new_DFT_calls'] == report['scheduler_mutations'] == 0
    assert report['maximum_new_interior_SCFs_per_chain'] == 7*(report['maximum_updates']+1)
    assert [r['target_strain'] for r in report['seeds']] == [.005, .01]
    for row in report['seeds']:
        images = read(output/row['label']/'seed.traj', ':')
        assert len(images) == 9 and all(a.calc is None for a in images)
        assert row['source_job'] == '28949435' and row['source_step'] == 9
        assert row['physical_input_sha256'] == CONTRACT
        assert row['retained_pattern_constraints_in_NEB'] == []
        assert not row['independent_forecast_test']
        validate_periodic_path_lift(images)
        for a in images:
            assert a.get_chemical_symbols() == ['Hf']*4+['O']*8
            np.testing.assert_allclose(a.cell.array[:2], images[0].cell.array[:2], atol=1e-12)
        params = json.loads((output/row['label']/'factory_parameters.json').read_bytes())
        assert set(params) == {'source_directory', 'seed_cache_records'}
        assert params['seed_cache_records'][1:8] == [None]*7
        old = (CASE.with_name('midpoint_reverse_E269_20261018')/'seeds/PO_flip_T_pattern_preserving/factory_parameters.json'
               if row['target_strain'] == .005 else
               CASE.with_name('clamped_plus1_preserving_terminal_E145_20261012')/'completed_HF/vcneb_summary.json')
        old_params = json.loads(old.read_bytes())
        if row['target_strain'] == .01: old_params = old_params['calculator_parameters']
        assert params['seed_cache_records'][0] == old_params['seed_cache_records'][0]
        assert params['seed_cache_records'][8] == old_params['seed_cache_records'][8]


def test_preparation_refuses_existing_namespace(prepared):
    output, _ = prepared
    with pytest.raises(FileExistsError): load('prepare').prepare(output)


def test_exact_held_resources_and_identity():
    module = load('control')
    good = ('JobId=123 JobState=PENDING Reason=JobHeldUser Partition=hfacnormal01 NumCPUs=32 '
            'NumTasks=1 CPUs/Task=32 NumNodes=1-1 TimeLimit=18:00:00 JobName=hfo2-E275-strain-0')
    module.validate_held(good, 0, '123')
    for raw in (good.replace('NumCPUs=32', 'NumCPUs=1'), good.replace('NumNodes=1-1', 'NumNodes=2'),
                good.replace('18:00:00', '1-00:00:00'), good.replace('JobId=123', 'JobId=124')):
        with pytest.raises(ValueError): module.validate_held(raw, 0, '123')


def test_unchanged_ordinary_runtime_and_native_mpirun():
    text = (CASE/'paths.slurm').read_text()
    for token in ('--steps 60', '--fmax .10', '--no-climb', 'mpirun -np 32',
                  '--time=18:00:00', '--cpus-per-task=32', 'I_MPI_HYDRA_BOOTSTRAP=fork',
                  'make_clamped_seed_cached_factory', '--require-continuous-periodic-lift'):
        assert token in text
