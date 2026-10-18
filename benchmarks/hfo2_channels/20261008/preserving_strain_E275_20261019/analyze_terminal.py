"""Same-strain native comparison; no DFT, model fitting or old-label replacement."""
from __future__ import annotations
import argparse
import csv
import importlib.util
import json
from pathlib import Path
import sys
import tarfile
import tempfile

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read
from examples.hfo2_fixed_input_factory import CONTRACT, read_fixed_hfo2_stru, same_ordered_geometry
from scripts.audit_hfo2_static_replica import audited_results, sha256
from scripts.analyze_hfo2_physical_holdout import analyze as original_data
from scripts.plot_hfo2_preserving_connection import restore as restore_zero, native_images, band as describe_chain
from scripts.plot_hfo2_mechanism import original_chart_images, path_coordinates, rotated_t_triplet

CASE = Path(__file__).resolve().parent
BASE = CASE.parent
SCALE = 5.12968067458423


def helper(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


def observation_images(folder):
    obs = json.loads((folder/'observation.json').read_bytes())
    if (obs['n_total_images'] != 9 or obs['formula_units'] != 4 or obs['pressure_GPa'] != 0
            or obs['climb'] or obs['fmax_target_eV_A'] != .10 or obs['cell_scale_A'] != SCALE
            or sha256(folder/'runtime_preflight.json') != obs['runtime_preflight_sha256']):
        raise ValueError('same ordinary nine-image P=0 fixed-physics observation required')
    images = []
    for i, row in enumerate(obs['raw_image_evaluations']):
        raw = folder/'raw'/f'image_{i:04d}'; poscar = folder/f'POSCAR_{i:02d}'
        if (row['image_index'] != i or sha256(poscar) != row['snapshot_POSCAR_sha256']
                or any(row['input_sha256'][n] != h for n, h in CONTRACT.items())
                or any(sha256(raw/n) != row['input_sha256'][n] for n in ('INPUT', 'KPT', 'STRU'))
                or sha256(raw/'OUT.ABACUS/running_scf.log') != row['raw_log_sha256']):
            raise ValueError('portable native input, ordered snapshot or log changed')
        image = read(poscar); parsed = audited_results(raw)
        expected = dict(energy=row['energy_eV_cell'], forces=row['forces_eV_A'], stress=row['stress_ASE_voigt_eV_A3'])
        if (not same_ordered_geometry(image, read_fixed_hfo2_stru(raw/'STRU'))
                or any(not np.allclose(parsed[k], expected[k], atol=1e-12, rtol=0) for k in parsed)):
            raise ValueError('native ordered E/F/stress does not replay')
        image.calc = SinglePointCalculator(image, **parsed); images.append(image)
    if len(images) != 9: raise ValueError('all nine native images required')
    return images, obs


def describe(images, obs, strain, kind, source_sha):
    t = read(BASE/'reference_variants/T.vasp')
    frame = json.loads((BASE/'clamped_endpoint_seeds/clamped_seed_manifest.json').read_bytes())['orientation']
    parent, patterns, _ = rotated_t_triplet(t)
    b = describe_chain(images, obs['mechanical_boundary']['reference_cell_A'], t, parent, patterns, frame, kind)
    _, _, chart, _ = path_coordinates(original_chart_images(images, frame, t), parent, patterns)
    if not np.isclose(b['fmax_eV_A'], obs['replayed_fmax_eV_A'], atol=1e-9, rtol=0):
        raise ValueError('full native band force differs from the exported observation')
    expected_plane = np.array(json.loads((BASE/'clamped_endpoint_seeds/clamped_seed_manifest.json').read_bytes())
                              ['unstrained_T_reference_cell_A'])[:2]*(1+strain)
    if not all(np.allclose(a.cell.array[:2], expected_plane, atol=1e-9, rtol=0) for a in images):
        raise ValueError('actual mechanical boundary does not match the declared strain')
    for i, row in enumerate(b['rows']):
        row.update(strain=strain, kind=kind, open_y_stretch_percent=float(100*(chart['deformation_gradients'][i, 1, 1]-1)))
    b.update(strain=strain, kind=kind, source_job=obs['source_job_id'], observation_sha256=source_sha,
             endpoint_native_records=[obs['raw_image_evaluations'][i] for i in (0, 8)])
    return b


def require_same_well(a, b, *, final=False):
    indices = (0, 1) if final else (0,)
    for i in indices:
        x, y = a['endpoint_native_records'][i], b['endpoint_native_records'][i]
        if (x['input_sha256'] != y['input_sha256'] or x['raw_log_sha256'] != y['raw_log_sha256']
                or x['energy_eV_cell'] != y['energy_eV_cell']):
            raise ValueError('same-strain competition requires the exact common native well/cache')


def reference_data():
    """Execute the comparison definitions on existing actual bands, not mock targets."""
    original = original_data(REPO); bands = []
    for row in original['bands']:
        if row['channel'] not in ('preserving', 'T'): continue
        folder = (REPO/row['rows'][0]['source_POSCAR']).parent
        images, obs = observation_images(folder)
        bands.append(describe(images, obs, row['strain'], 'registered_preserving' if row['channel'] == 'preserving'
                              else 'T_escape', sha256(folder/'observation.json')))
    with tempfile.TemporaryDirectory(prefix='E275-zero-') as td:
        folder, _ = restore_zero(9, Path(td), terminal=True)
        images, obs, _ = native_images(folder, 9)
        bands.append(describe(images, obs, 0., 'released_seed', sha256(folder/'observation/observation.json')))
    previous = helper('E275_midpoint_archive', BASE/'midpoint_reverse_E269_20261018/analyze_completed.py')
    with tempfile.TemporaryDirectory(prefix='E275-midpoint-') as td:
        exported, _ = previous.restore(Path(td))
        folder = exported/'preserving/step_0026'
        images, obs = observation_images(folder)
        bands.append(describe(images, obs, .005, 'upper_seed_midpoint', sha256(folder/'observation.json')))
    for strain in (0., .005, .01):
        same = [b for b in bands if b['strain'] == strain]
        p = next(b for b in same if b['kind'] == 'registered_preserving')
        for b in same: require_same_well(p, b, final=b['kind'] != 'T_escape')
    return dict(bands=bands, original_forecast_manifest_sha256=original['forecast_manifest_sha256'],
                old_labels_and_scores_unchanged=True, new_DFT_calls=0)


def read_terminal(label, destination):
    source = CASE/'completed'/label
    outer = json.loads((source/'collection_receipt.json').read_bytes())
    archive = source/'native_terminal.tar.gz'
    collector = helper('E275_offline_admission', CASE/'collect_terminal.py')
    if (outer['job_id'] != collector.JOBS[label] or outer['label'] != label
            or outer['bundle_pin'] != collector.PIN or sha256(archive) != outer['archive_sha256']
            or outer['new_DFT_calls'] != 0 or outer['scheduler_mutations'] != 0
            or not outer['all_runtime_calls_exported']):
        raise ValueError('complete original native archive/identity required')
    with tarfile.open(archive) as tar:
        for member in tar.getmembers():
            p = Path(member.name)
            if (not p.parts or p.parts[0] != 'native_terminal' or '..' in p.parts
                    or not(member.isfile() or member.isdir())
                    or not(destination/p).resolve().is_relative_to(destination.resolve())):
                raise ValueError('unsafe native archive member')
        tar.extractall(destination)
    folder = destination/'native_terminal'
    files = json.loads((folder/'files.json').read_bytes())
    if len(files) != outer['exported_files']: raise ValueError('complete archive payload count differs')
    for name, digest in files.items():
        p = folder/name
        if not p.resolve().is_relative_to(folder.resolve()) or sha256(p) != digest:
            raise ValueError('archived native bytes changed')
    receipt = json.loads((folder/'receipt.json').read_bytes())
    summary = json.loads((folder/'vcneb_summary.json').read_bytes())
    preflight = json.loads((folder/'vcneb_preflight.json').read_bytes())
    row = receipt['scheduler_row_at_export']; job = receipt['job_id']; step = receipt['step']
    collector.terminal_row('|'.join(row), job)
    expected = collector.admit(job, row, summary, step, False, preflight, receipt['producer_sha256'])
    if (job != outer['job_id'] or receipt['bundle_pin'] != collector.PIN
            or receipt['producer_sha256'] != sha256(CASE/'paths.slurm')
            or receipt['new_SCF_calls'] != expected or len(receipt['calls']) != expected
            or receipt['ordinary_converged'] != outer['ordinary_converged']):
        raise ValueError('terminal metadata/call cohort differs')
    counts = {f'image_{i:04d}': 0 for i in range(1, 8)}; seconds = []
    for record in receipt['calls']:
        relative = Path(record['relative_directory']); raw = folder/'all_calls'/relative
        if relative.parts[0] not in counts: raise ValueError('fixed endpoint was recalculated')
        counts[relative.parts[0]] += 1
        call = json.loads((raw/'call_audit.json').read_bytes())
        if (sha256(raw/'call_audit.json') != record['call_audit_sha256']
                or call['input_sha256'] != record['input_sha256']
                or any(call['input_sha256'][n] != h for n, h in CONTRACT.items())
                or any(sha256(raw/n) != call['input_sha256'][n] for n in ('INPUT', 'KPT', 'STRU'))
                or sha256(raw/'OUT.ABACUS/running_scf.log') != record['raw_log_sha256']):
            raise ValueError('native call input/log changed')
        parsed = audited_results(raw)
        if any(not np.allclose(parsed[k], call['results'][k], atol=1e-12, rtol=0) for k in parsed):
            raise ValueError('native intermediate E/F/stress differs')
        if call['elapsed_seconds'] != record['seconds'] or not np.isfinite(record['seconds']) or record['seconds'] <= 0:
            raise ValueError('positive unchanged timing required')
        seconds.append(record['seconds'])
    if (any(n != step+1 for n in counts.values())
            or not np.isclose(sum(seconds)*32/3600, receipt['native_proxy_core_hours'], atol=1e-10, rtol=0)
            or not np.isclose(float(row[7])*32/3600, receipt['allocated_core_hours'], atol=1e-10, rtol=0)):
        raise ValueError('complete per-image cohort or cost differs')
    images, obs = observation_images(folder/'observation')
    band = describe(images, obs, outer['target_strain'], 'released_seed', sha256(folder/'observation/observation.json'))
    if (band['ordinary_residual_pass'] != receipt['ordinary_converged']
            or not np.isclose(band['sampled_polygon_peak_meV_fu'], receipt['sampled_peak_meV_fu'], atol=1e-9, rtol=0)):
        raise ValueError('native terminal peak/residual does not replay')
    return band, dict(job_id=job, calls=expected, step=step,
        ordinary_converged=receipt['ordinary_converged'], archive_sha256=outer['archive_sha256'],
        native_proxy_core_hours=receipt['native_proxy_core_hours'], allocated_core_hours=receipt['allocated_core_hours'],
        all_new_native_EFS_reparsed_locally=True, full_orbital_bytes_rechecked_locally=False)


def summarize(bands):
    rows = []
    for strain in (0., .005, .01):
        group = {b['kind']: b for b in bands if b['strain'] == strain}
        old, new, escape = (group[k] for k in ('registered_preserving', 'released_seed', 'T_escape'))
        require_same_well(old, new, final=True); require_same_well(new, escape)
        rows.append(dict(strain=strain, original_preserving_meV_fu=old['sampled_polygon_peak_meV_fu'],
            released_seed_meV_fu=new['sampled_polygon_peak_meV_fu'], T_meV_fu=escape['sampled_polygon_peak_meV_fu'],
            old_T_minus_preserving_meV_fu=escape['sampled_polygon_peak_meV_fu']-old['sampled_polygon_peak_meV_fu'],
            new_T_minus_preserving_meV_fu=escape['sampled_polygon_peak_meV_fu']-new['sampled_polygon_peak_meV_fu'],
            new_fmax_eV_A=new['fmax_eV_A'], new_ordinary_residual_pass=new['ordinary_residual_pass'],
            new_peak_image=new['sampled_peak_image'], new_peak_polar_y_A=new['rows'][new['sampled_peak_image']]['polar_proxy_y_A'],
            new_peak_open_y_percent=new['rows'][new['sampled_peak_image']]['open_y_stretch_percent']))
    return rows


def analyze(output):
    output = Path(output)
    if output.exists(): raise FileExistsError(output)
    # Missing/live jobs remain missing: no partial target, selection or relaunch.
    for label in ('strain_0050', 'strain_0100'):
        if not (CASE/'completed'/label/'collection_receipt.json').is_file():
            raise FileNotFoundError(f'await original terminal archive: {label}')
    data = reference_data(); costs = []
    for label in ('strain_0050', 'strain_0100'):
        with tempfile.TemporaryDirectory(prefix='E275-terminal-') as td:
            b, cost = read_terminal(label, Path(td)); data['bands'].append(b); costs.append(cost)
    data.update(outcomes=summarize(data['bands']), terminal_cohorts=costs,
        status='finite_same_strain_connection_comparison', physical_error_bound=None,
        global_MEP_or_TS_certified=False, independent_prediction_test=False)
    output.mkdir(parents=True)
    with (output/'analysis.json').open('x', encoding='utf-8', newline='\n') as h:
        json.dump(data, h, indent=2, allow_nan=False); h.write('\n')
    for name, rows in (('outcomes.csv', data['outcomes']),
                       ('profiles.csv', [r for b in data['bands'] for r in b['rows']])):
        with (output/name).open('x', encoding='utf-8', newline='') as h:
            writer = csv.DictWriter(h, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    return data


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    r = analyze(p.parse_args().output); print(json.dumps(r['outcomes']))
