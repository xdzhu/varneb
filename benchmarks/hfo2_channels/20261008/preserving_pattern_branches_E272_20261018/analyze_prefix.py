"""Native replay and full-polygon connection check; no calculator calls.

The selected J geometry remains an ordinary, free interior image. Its local
energy lowering is not subtracted from a path barrier or an old NEB residual.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import sys
import tarfile
import tempfile

REPOSITORY = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPOSITORY))
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read, write
from examples.hfo2_fixed_input_factory import CONTRACT, read_fixed_hfo2_stru, same_ordered_geometry
from scripts.audit_hfo2_static_replica import audited_results, sha256
from vcneb import VCNEB, clamped_plane_vcneb_boundary, validate_path_geometry, validate_periodic_path_lift

CASE = Path(__file__).resolve().parent
BASE = CASE.parent
ARCHIVE_SHA = '3570708a58716934591b13febc82ec6e42ef647db13ae0b1b3f78338af6e9501'
BUNDLE_PIN = 'd12af56eac84edf28e536c69b524d12681f876db1918c889fbf4d89eb59e06c2'
OBSERVATION_SHA = 'babd26a259339ac4a8fb9f0154895b7292294e97bdc8f2ad478fbfc2312e8995'
OBSERVATION = BASE/'clamped_held_flips_E100_20261011/zero_strain_preserving/observations/step_0015'
SCALE = 5.12968067458423


def restore(destination):
    archive = CASE/'prefix_0_minus.tar.gz'
    receipt = json.loads((CASE/'prefix_0_minus_receipt.json').read_bytes())
    if sha256(archive) != ARCHIVE_SHA or receipt['archive_sha256'] != ARCHIVE_SHA:
        raise ValueError('selected native prefix archive changed')
    if (receipt['selected_native_SCFs'] != 4 or receipt['whole_job_completed']
            or receipt['all_attempts_exported'] or not receipt['partial_signed_branch_only']):
        raise ValueError('partial-prefix scope differs')
    with tarfile.open(archive) as tar:
        for member in tar.getmembers():
            name = Path(member.name)
            if (name.is_absolute() or '..' in name.parts or not name.parts
                    or name.parts[0] != 'prefix_0_minus' or not (member.isfile() or member.isdir())
                    or not (destination/name).resolve().is_relative_to(destination.resolve())):
                raise ValueError('unsafe selected native archive member')
        tar.extractall(destination)
    root = destination/'prefix_0_minus'
    files = json.loads((root/'files.json').read_bytes())
    if len(files) != receipt['exported_files']:
        raise ValueError('native prefix manifest count differs')
    for name, digest in files.items():
        path = root/name
        if not path.resolve().is_relative_to(root.resolve()) or sha256(path) != digest:
            raise ValueError('selected native bytes changed')
    return root


def setup(index=0):
    if index not in (0, 2):
        raise ValueError('only original zero and +1% preserving neighborhoods admitted')
    manifest = CASE/'executed_bundle_files.json'
    if sha256(manifest) != BUNDLE_PIN:
        raise ValueError('executed E272 physical protocol changed')
    for name, digest in json.loads(manifest.read_bytes()).items():
        if sha256(CASE/name) != digest:
            raise ValueError('executed E272 source/input changed')
    original = BASE/'matched_pattern_release_E270_20261018'
    files = json.loads((original/'executed_bundle_files.json').read_bytes())
    if sha256(original/'study.py') != files['study.py']:
        raise ValueError('original chart/native replay source changed')
    spec = importlib.util.spec_from_file_location('E272_prefix_original_chart', original/'study.py')
    study = importlib.util.module_from_spec(spec); spec.loader.exec_module(study)
    seed = CASE/'prepared/seeds'/str(index)
    data = json.loads((seed/'seed.json').read_bytes())
    expected = (0., '28727382') if index == 0 else (.01, '28812952')
    if (data['channel'] != 'preserving' or (data['substrate_strain'],data['source_job']) != expected
            or index == 0 and data['observation_sha256'] != OBSERVATION_SHA):
        raise ValueError('original registered preserving peak required')
    chart = study.chart_for(seed)
    e0, g0 = study.cached_frozen(seed, chart)
    return chart, data, e0, g0


def replay_native(root, chart, coordinates):
    x = np.asarray(coordinates)
    if x.shape != (39,) or not np.isfinite(x).all():
        raise ValueError('finite original-chart coordinate required')
    matches = [p for p in root.glob('evaluations/eval-*/result.json')
               if json.loads(p.read_bytes())['coordinates'] == coordinates]
    if len(matches) != 1:
        raise ValueError('one exact native result is required for each coordinate')
    path = matches[0]; native = json.loads(path.read_bytes()); v = native['validation']
    raw = path.parent/'scf_000000'
    if (native['status'] != 'complete' or v['genuine_MPI_ranks'] != 32
            or any(v['input_sha256'][n] != h for n, h in CONTRACT.items())
            or any(sha256(raw/n) != v['input_sha256'][n] for n in ('INPUT', 'KPT', 'STRU'))
            or sha256(raw/'OUT.ABACUS/running_scf.log') != v['raw_log_sha256']
            or sha256(raw/'call_audit.json') != v['call_audit_sha256']):
        raise ValueError('native SCF input, MPI or output identity differs')
    atoms = chart.to_atoms(x)
    if not same_ordered_geometry(atoms, read_fixed_hfo2_stru(raw/'STRU')):
        raise ValueError('coordinate does not match its actual ordered native input')
    results = audited_results(raw)
    audit = json.loads((raw/'call_audit.json').read_bytes())
    if any(not np.allclose(results[k], audit['results'][k], atol=1e-12, rtol=0) for k in results):
        raise ValueError('native energy/forces/stress differ from call record')
    atoms.calc = SinglePointCalculator(atoms, **results)
    energy, gradient = chart.energy_and_gradient(x, atoms)
    if (energy != native['enthalpy_eV']
            or not np.allclose(gradient, native['gradient'], atol=2e-12, rtol=0)):
        raise ValueError('independent physical-gradient replay differs')
    return atoms, energy, gradient, dict(native_result_sha256=sha256(path),
        native_record_relative_path=path.relative_to(root).as_posix(),
        native_log_sha256=v['raw_log_sha256'], native_input_sha256=v['input_sha256'],
        full_orbital_pseudopotential_bytes_rechecked_locally=False)


def review_branch(root, chart, data, e0):
    summary = json.loads((root/'treatments.json').read_bytes())
    if (summary['index'] != 0 or summary['strain'] != 0 or summary['sign'] != 'minus'
            or summary['total_new_attempts'] != 15 or summary['frozen_energy_eV_cell'] != e0
            or summary['stationary_or_stability_certified'] or not summary['old_path_labels_unchanged']):
        raise ValueError('the completed signed zero-strain branch is required')
    q = np.asarray(data['orthonormal_39_by_4']); c = np.asarray(data['covectors_39_by_4'])
    translation = np.zeros((39, 3)); translation[:36] = np.tile(np.eye(3), (12, 1))/np.sqrt(12)
    rows = {}; candidates = {}
    for arm, attempts in (('A', 12), ('J', 1)):
        treatment = json.loads((root/f'treatment_{arm}.json').read_bytes())
        if (treatment != summary['treatments'][0 if arm == 'A' else 1]
                or treatment['treatment'] != arm or treatment['new_attempts'] != attempts
                or treatment['status'] != 'gradient_candidate_stability_unchecked'):
            raise ValueError('declared arm identity, result or finite budget differs')
        candidate = treatment['candidate']; x = np.asarray(candidate['coordinates'])
        if arm == 'A' and np.any(x[36:]):
            raise ValueError('A must retain an exactly fixed cell')
        atoms, energy, g, evidence = replay_native(root, chart, candidate['coordinates'])
        n = 36 if arm == 'A' else 39
        released = g[:n]-q[:n]@(q[:n].T@g[:n])-translation[:n]@(translation[:n].T@g[:n])
        residual = float(np.linalg.norm(released)); defect = float(np.max(np.abs(c.T@x)))
        if (energy != candidate['energy_eV_cell'] or residual > .03 or defect > 1e-10
                or not np.isclose(residual, candidate['released_gradient_norm_eV_A'], atol=2e-12, rtol=0)
                or candidate['stationary_TS_or_stability_certified']):
            raise ValueError('local candidate target, energy or physical constraints differ')
        rows[arm] = dict(energy_eV_cell=energy, dE_minus_F_meV_fu=250*(energy-e0),
            released_gradient_norm_eV_A=residual, full_allowed_gradient_norm_eV_A=float(np.linalg.norm(g)),
            retained_pattern_defect_A=defect, **evidence)
        candidates[arm] = atoms
    check = json.loads((root/'known_direction_check.json').read_bytes())
    probes = [replay_native(root, chart, x) for x in check['coordinates']]
    direction = (np.asarray(check['coordinates'][1])-np.asarray(check['coordinates'][0]))/.02
    prep = json.loads((CASE/'prepared/seeds/0/branch_preparation.json').read_bytes())
    if (check != summary['known_direction_check'] or check['new_attempts'] != 2 or check['step_A'] != .01
            or check['whole_stability_certified']
            or not np.allclose(direction, np.asarray(prep['signed_starts'][1])/.1, atol=2e-12, rtol=0)):
        raise ValueError('only the registered direction at A may be checked')
    kg = float((probes[1][2]-probes[0][2])@direction/.02)
    ke = float((probes[1][1]-2*rows['A']['energy_eV_cell']+probes[0][1])/.01**2)
    if (not np.isclose(kg, check['gradient_secant_eV_A2'], atol=2e-12, rtol=0)
            or not np.isclose(ke, check['energy_secant_eV_A2'], atol=2e-12, rtol=0)):
        raise ValueError('native directional secant differs')
    return candidates['J'], dict(arms=rows, directional_secants_at_A_eV_A2=dict(gradient=kg, energy=ke),
        full_stability_at_A_or_J_certified=False, local_native_SCFs_reparsed=4)


def original_images():
    path = OBSERVATION/'observation.json'
    if sha256(path) != OBSERVATION_SHA:
        raise ValueError('original preserving nine-image observation changed')
    obs = json.loads(path.read_bytes())
    if (obs['source_job_id'] != '28727382' or obs['n_total_images'] != 9 or obs['formula_units'] != 4
            or obs['pressure_GPa'] != 0 or obs['climb'] or obs['cell_scale_A'] != SCALE
            or obs['k_eV_A2'] != .2 or obs['fmax_target_eV_A'] != .10
            or sha256(OBSERVATION/'runtime_preflight.json') != obs['runtime_preflight_sha256']):
        raise ValueError('ordinary same-boundary input/driver contract differs')
    images = []
    for i, record in enumerate(obs['raw_image_evaluations']):
        poscar = OBSERVATION/f'POSCAR_{i:02d}'; raw = OBSERVATION/'raw'/f'image_{i:04d}'
        if (record['image_index'] != i or sha256(poscar) != record['snapshot_POSCAR_sha256']
                or any(record['input_sha256'][n] != h for n, h in CONTRACT.items())
                or any(sha256(raw/n) != record['input_sha256'][n] for n in ('INPUT', 'KPT', 'STRU'))
                or sha256(raw/'OUT.ABACUS/running_scf.log') != record['raw_log_sha256']):
            raise ValueError('original image native input/output bytes differ')
        image = read(poscar, format='vasp'); results = audited_results(raw)
        expected = dict(energy=record['energy_eV_cell'], forces=record['forces_eV_A'],
                        stress=record['stress_ASE_voigt_eV_A3'])
        if (not same_ordered_geometry(image, read_fixed_hfo2_stru(raw/'STRU'))
                or any(not np.allclose(results[k], expected[k], atol=1e-12, rtol=0) for k in results)):
            raise ValueError('ordered image or independent native E/F/stress replay differs')
        image.calc = SinglePointCalculator(image, **results); images.append(image)
    return images, obs


def fixed_lift(original, native_frozen, candidate):
    """Transport one fixed integer lift, with no new atom mapping or MIC choice."""
    if not same_ordered_geometry(original, native_frozen):
        raise ValueError('native F differs from the original ordered peak')
    delta = original.get_scaled_positions(wrap=False)-native_frozen.get_scaled_positions(wrap=False)
    integers = np.rint(delta)
    if not np.allclose(delta, integers, atol=1e-10, rtol=0):
        raise ValueError('original/native F gauges are not a fixed integer lift')
    lifted = candidate.copy()
    lifted.set_scaled_positions(candidate.get_scaled_positions(wrap=False)+integers)
    if not same_ordered_geometry(lifted, candidate):
        raise ValueError('integer lift changed the actual ordered periodic geometry')
    lifted.calc = SinglePointCalculator(lifted, **candidate.calc.results)
    return lifted, integers.astype(int).tolist()


def polygon(images, reference_cell):
    boundary = clamped_plane_vcneb_boundary(12, reference_cell, allow_tilt=True)
    boundary.validate_images(images); validate_periodic_path_lift(images)
    geometry = validate_path_geometry(images, cell_scale=SCALE, minimum_distance=1.6, maximum_deformation=.25)
    chain = VCNEB(images, cell_scale=SCALE, k=.2, climb=False, pressure=0., **boundary.vcneb_kwargs(images))
    forces = chain.get_forces().reshape(7, -1, 3)
    if forces.shape != (7, 15, 3):
        raise ValueError('registered atomic plus cell-row force packing differs')
    per_image = np.linalg.norm(forces, axis=2).max(axis=1)
    fmax = float(per_image.max()); relative = (chain.enthalpies-chain.enthalpies[0])*250
    return dict(replayed_fmax_eV_A=fmax, ordinary_residual_pass=bool(fmax <= .10),
        per_interior_image_fmax_eV_A=per_image.tolist(),
        per_interior_image_atomic_fmax_eV_A=np.linalg.norm(forces[:, :12], axis=2).max(axis=1).tolist(),
        per_interior_image_cell_row_fmax_eV_A=np.linalg.norm(forces[:, 12:], axis=2).max(axis=1).tolist(),
        energy_eV_cell=chain.enthalpies.tolist(), relative_energy_meV_fu=relative.tolist(),
        sampled_peak_image=int(np.argmax(relative)), sampled_polygon_peak_meV_fu=float(relative.max()),
        interior_local_peak_images=[i for i in range(1, 8)
                                   if relative[i] > max(relative[i-1], relative[i+1])],
        extended_arc_A=chain.reaction_coordinate().tolist(), geometry_gate=geometry)


def analyze():
    chart, data, e0, _ = setup()
    with tempfile.TemporaryDirectory(prefix='varneb_P_native_prefix_') as temporary:
        root = restore(Path(temporary)); candidate, branch = review_branch(root, chart, data, e0)
    images, obs = original_images(); reference = obs['mechanical_boundary']['reference_cell_A']
    original = polygon(images, reference)
    if (abs(original['replayed_fmax_eV_A']-obs['replayed_fmax_eV_A']) > 7e-7
            or not np.allclose(original['relative_energy_meV_fu'], obs['relative_enthalpy_meV_fu'], atol=1e-8, rtol=0)):
        raise ValueError('original same-boundary chain replay differs')
    replacement, lift = fixed_lift(images[4], chart.reference_atoms, candidate)
    images[4] = replacement
    changed = polygon(images, reference)
    report = dict(status='selected_zero_minus_native_branch_and_full_polygon_replayed',
        original=original, replaced_central_image=changed, signed_branch=branch,
        replaced_image=4, fixed_per_atom_integer_lift=lift, new_DFT_calls=0,
        n_total_images=9, n_fixed_endpoints=2, n_free_interior_images=7, formula_units=4,
        pressure_GPa=0., climb=False, k_eV_A2=.2, cell_scale_A=SCALE, fmax_target_eV_A=.10,
        archive_sha256=ARCHIVE_SHA, original_observation_sha256=OBSERVATION_SHA,
        original_labels_and_forecast_scores_unchanged=True, whole_E272_job_completed_here=False,
        interpolating_segments_evaluated=False, stationary_TS_or_MEP_or_corrected_barrier_certified=False,
        interpretation='single-image replacement is a connection diagnostic, not a relaxed new path',
        limitations=['Only four selected branch SCFs, not all 15 attempts, are in this prefix.',
                     'The positive secants are one direction at A, not a full A/J stability test.',
                     'Every polygon tangent, spring and allowed atomic/cell force is recomputed.',
                     'A local lowering cannot be subtracted from an old path barrier.',
                     'J remains a free interior image; endpoints and physical settings are unchanged.'])
    return report, images


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); report, images = analyze(); args.output.mkdir(parents=True, exist_ok=False)
    with (args.output/'analysis.json').open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(report, handle, indent=2, allow_nan=False); handle.write('\n')
    write(args.output/'evaluated_polygon.traj', images)
    for index, image in enumerate(images):
        write(args.output/f'POSCAR_{index:02d}', image, format='vasp', direct=True, vasp5=True)
    print(json.dumps({key: report[key] for key in ('status', 'new_DFT_calls', 'original', 'replaced_central_image')}))
