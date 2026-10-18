"""Upper-to-midpoint geometry and bounded dispatch; no DFT or scheduler."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
from ase import Atoms
from ase.io import read

from scripts.prepare_hfo2_midpoint_reverse import CASE, OLD, prepare, transfer_endpoint_lift
from examples.hfo2_fixed_input_factory import same_ordered_geometry
from vcneb import validate_periodic_path_lift

ROOT = Path(__file__).resolve().parents[1]


def test_endpoint_lift_keeps_integer_gauge_without_permuting():
    anchor = Atoms(['Hf', 'O'], cell=np.eye(3)*8,
                   scaled_positions=[[.1, .2, .3], [.4, .5, .6]], pbc=True)
    source = anchor.copy()
    source.set_scaled_positions(anchor.get_scaled_positions(wrap=False)+[[1.001, 0, 0], [0, -.999, 0]])
    target = anchor.copy(); target.positions[:, 2] += .003
    result, registration = transfer_endpoint_lift(source, target, anchor)
    np.testing.assert_array_equal(registration['integer_lift'], [[1, 0, 0], [0, -1, 0]])
    np.testing.assert_allclose(result.get_scaled_positions(wrap=False),
        target.get_scaled_positions(wrap=False)+[[1, 0, 0], [0, -1, 0]], atol=1e-12)
    assert same_ordered_geometry(result, target)
    assert result.calc is None


def test_ambiguous_parent_motion_or_species_is_not_remapped():
    anchor = Atoms(['Hf', 'O'], cell=np.eye(3)*8,
                   scaled_positions=[[.1, .2, .3], [.4, .5, .6]], pbc=True)
    wrong = anchor.copy(); wrong.positions[0, 0] += 2.4
    with pytest.raises(ValueError, match='registration'):
        transfer_endpoint_lift(wrong, anchor, anchor)
    wrong = anchor.copy(); wrong.set_chemical_symbols(['O', 'Hf'])
    with pytest.raises(ValueError, match='species'):
        transfer_endpoint_lift(wrong, anchor, anchor)


@pytest.fixture(scope='module')
def generated(tmp_path_factory):
    output = tmp_path_factory.mktemp('reverse')/'seeds'
    return output, prepare(ROOT, output)


def test_actual_two_sources_target_boundary_and_same_endpoint_caches(generated):
    output, report = generated
    assert len(report['seeds']) == 2 and report['new_DFT_calls'] == 0
    assert not report['independent_forecast_test']
    for row in report['seeds']:
        folder = output/row['channel']
        images = read(folder/'seed.traj', index=':')
        original = read(ROOT/OLD/'path_seeds'/row['channel']/'seed.traj', index=':')
        assert len(images) == 9 and all(a.calc is None for a in images)
        assert row['source_strain'] == .01 and row['target_strain'] == .005
        assert row['source_job'] == ('28812952' if 'preserving' in row['channel'] else '28709788')
        validate_periodic_path_lift(images)
        assert all(same_ordered_geometry(a, b) for a, b in zip((images[0], images[-1]), (original[0], original[-1])))
        for image in images:
            np.testing.assert_allclose(image.cell.array[:2], original[0].cell.array[:2], atol=1e-12)
        new = json.loads((folder/'factory_parameters.json').read_bytes())
        old = json.loads((ROOT/OLD/'path_seeds'/row['channel']/'factory_parameters.json').read_bytes())
        assert new == old
        assert new['seed_cache_records'][1:8] == [None]*7
        assert not same_ordered_geometry(images[4], original[4])


def test_existing_output_not_overwritten(generated):
    output, _ = generated
    original = (output/'preparation.json').read_bytes()
    with pytest.raises(FileExistsError):
        prepare(ROOT, output)
    assert (output/'preparation.json').read_bytes() == original


def test_exact_held_resources_and_handle_required():
    spec = importlib.util.spec_from_file_location('E269_control', ROOT/CASE/'control.py')
    control = importlib.util.module_from_spec(spec); spec.loader.exec_module(control)
    good = ('JobId=123 JobState=PENDING Reason=JobHeldUser Partition=hfacnormal01 '
            'NumCPUs=32 NumTasks=1 CPUs/Task=32 NumNodes=1-1 TimeLimit=1-00:00:00 '
            'JobName=hfo2-E269-midpoint-0')
    control.validate_held(good, 0, '123')
    for bad in (good.replace('NumCPUs=32', 'NumCPUs=1'), good.replace('NumNodes=1-1', 'NumNodes=2'),
                good.replace('JobId=123', 'JobId=124'), good.replace('hfacnormal01', 'other')):
        with pytest.raises(ValueError):
            control.validate_held(bad, 0, '123')


def test_native_slurm_recipe_stays_ordinary_and_uses_mpirun32():
    text = (ROOT/CASE/'paths.slurm').read_text()
    assert '--no-climb' in text and '--fmax .10' in text
    assert '--steps 80' in text and '--time=24:00:00' in text
    assert '--cpus-per-task=32' in text and 'mpirun -np 32' in text
    assert 'I_MPI_HYDRA_BOOTSTRAP=fork' in text


def test_readonly_inspector_cli_has_no_stdlib_shadow_or_cwd_dependency(tmp_path):
    helper = ROOT/CASE/'inspect_paths.py'
    assert not helper.with_name('inspect.py').exists()
    result = subprocess.run([sys.executable, str(helper), '--help'], cwd=tmp_path,
                            capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr
    assert '--audit-prefix' in result.stdout


def test_six_physical_neighborhoods_replay_allowed_not_raw_reaction_force():
    spec = importlib.util.spec_from_file_location('E269_neighborhoods', ROOT/CASE/'neighborhoods.py')
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    report = module.register(ROOT)
    rows = report['rows']
    assert len(rows) == 6 and report['new_DFT_calls'] == report['new_curvature_probes'] == 0
    assert {(r['channel'], r['substrate_strain']) for r in rows} == {
        (c, s) for c in ('preserving', 'T') for s in (0., .005, .01)}
    assert [r['sampled_peak_meV_fu'] for r in rows] == pytest.approx([
        96.8281830028, 98.5308769164, 125.4244311581, 84.4166238812, 90.4181783358, 96.7636503947])
    for row in rows:
        assert not row['stationary_bottleneck_certified'] and not row['Hessian_index_certified']
        assert row['peak_allowed_physical_norm_eV_A'] >= row['peak_allowed_physical_max_vector_eV_A']
        assert abs(row['peak_allowed_tangential_force_eV_A']) <= row['peak_allowed_physical_norm_eV_A']
        assert row['ordinary_chain_fmax_eV_A'] < .10
    assert abs(rows[1]['peak_allowed_tangential_force_eV_A']) < .002
    assert rows[1]['peak_allowed_physical_norm_eV_A'] > .2
