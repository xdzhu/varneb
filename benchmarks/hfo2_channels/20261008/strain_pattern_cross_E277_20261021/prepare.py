"""Prepare exactly two cross controls from existing native J candidates."""
from __future__ import annotations
import argparse
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import tempfile
import numpy as np
from ase.io import read, write

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))
from scripts.audit_hfo2_static_replica import sha256
from scripts.plot_hfo2_mechanism import original_chart_images, path_coordinates, rotated_t_triplet
from vcneb import validate_path_geometry

CASE = Path(__file__).resolve().parent
BASE = CASE.parent
OLD = BASE/'preserving_pattern_branches_E272_20261018'
NEW = BASE/'preserving_local_completion_E276_20261019'
sys.path.insert(0, str(CASE))
from policy import POLICY


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def save(path, value):
    with Path(path).open('x', encoding='utf-8', newline='\n') as f:
        json.dump(value, f, indent=2, allow_nan=False); f.write('\n')


def observables(atoms):
    t = read(BASE/'reference_variants/T.vasp')
    orientation = json.loads((BASE/'clamped_endpoint_seeds/clamped_seed_manifest.json').read_bytes())['orientation']
    parent, patterns, _ = rotated_t_triplet(t)
    images = original_chart_images([atoms, atoms], orientation, t)
    q, polar, _, _ = path_coordinates(images, parent, patterns)
    return np.r_[q[0, [0, 2, 1]], polar[0, 2]], polar[0], images[0].cell.array


def transported_start(target_chart, donor_chart, donor_atoms):
    """Register once from the original frozen pair, then transport the donor."""
    ref = target_chart.reference_atoms
    sf = ref.get_scaled_positions(wrap=False)
    sd = donor_chart.reference_atoms.get_scaled_positions(wrap=False)
    shift = np.rint(sf-sd)
    if np.max(np.abs(sd+shift-sf)) > .2:
        raise ValueError('frozen pair does not admit the declared short registered lift')
    start = donor_atoms.copy(); start.calc = None
    cell = start.cell.array.copy(); cell[:2] = ref.cell.array[:2]
    start.set_cell(cell, scale_atoms=True)
    s = donor_atoms.get_scaled_positions(wrap=False)+shift
    translation = np.mean(s-sf, axis=0)
    start.set_scaled_positions(s-translation)
    x = target_chart.from_atoms(start)
    if np.linalg.norm(x) >= 1.:
        raise ValueError('cross start lies outside the predeclared numerical window')
    report = validate_path_geometry([ref, start], cell_scale=5.12968067458423,
        minimum_distance=1.6, maximum_deformation=.25)
    if not report['valid']:
        raise ValueError('cross start geometry rejected')
    return start, x, dict(integer_shift_by_atom=shift.astype(int).tolist(),
        removed_fractional_translation=translation.tolist(), geometry=report)


def prepare(output):
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(output)
    output.mkdir(parents=True, exist_ok=True)
    reader = load(OLD/'analyze_prefix.py', 'E277_native_reader')
    old_archive = load(OLD/'analyze_terminal.py', 'E277_old_archive')
    new_archive = load(NEW/'analyze_terminal.py', 'E277_new_archive')
    diagonal = {}; charts = {}; metadata = {}; geometries = {}
    for label, index in (('q0', 0), ('q1', 2)):
        chart, data, _, _ = reader.setup(index)
        charts[label], metadata[label] = chart, data
        with tempfile.TemporaryDirectory(prefix='E277-diagonal-') as td:
            if label == 'q0':
                root, receipt = old_archive.restore(Path(td)); folder = root/'0/minus'
            else:
                root, receipt = new_archive.restore('minus', Path(td)); folder = root
            candidate = json.loads((folder/'treatment_J.json').read_bytes())['candidate']
            atoms, energy, g, evidence = reader.replay_native(folder, chart, candidate['coordinates'])
            q = np.asarray(data['orthonormal_39_by_4'])
            tr = chart.rigid_translation_directions()/np.sqrt(12)
            residual = float(np.linalg.norm(g-q@(q.T@g)-tr@(tr.T@g)))
            if residual > .03 or energy != candidate['energy_eV_cell']:
                raise ValueError('existing qualified joint candidate required')
            qvalues, polar, _ = observables(atoms)
            fvalues, _, _ = observables(chart.reference_atoms)
            if not np.allclose(qvalues, fvalues, atol=1e-12, rtol=0):
                raise ValueError('diagonal candidate changed its four physical observables')
            diagonal[label] = dict(strain=data['substrate_strain'], energy_eV_cell=energy,
                coordinates=candidate['coordinates'], held_values_A=qvalues.tolist(),
                polar_proxy_A=polar.tolist(), released_gradient_eV_A=residual,
                archive_sha256=receipt['archive_sha256'], **evidence)
            geometries[label] = atoms
            native = folder/evidence['native_record_relative_path']
            shutil.copytree(native.parent, output/'diagonals'/label/'native')
            write(output/'diagonals'/label/'candidate.traj', atoms)
    rows = []
    for name, target, donor in (('zero_q1', 'q0', 'q1'), ('plus1_q0', 'q1', 'q0')):
        chart, data = charts[target], metadata[target]
        start, x, registration = transported_start(chart, charts[donor], geometries[donor])
        held, _, _ = observables(start)
        c = np.asarray(data['covectors_39_by_4'])
        expected = np.asarray(diagonal[donor]['held_values_A'])-diagonal[target]['held_values_A']
        if (not np.allclose(held, diagonal[donor]['held_values_A'], atol=1e-12, rtol=0)
                or not np.allclose(c.T@x, expected, atol=1e-12, rtol=0)):
            raise ValueError('transport changed the donor physical-observable vector')
        folder = output/name
        shutil.copytree(OLD/'prepared/seeds'/str(0 if target == 'q0' else 2), folder/'seed')
        prep = json.loads((folder/'seed/branch_preparation.json').read_bytes())
        row = dict(name=name, target_pattern=target, held_pattern=donor,
            substrate_strain=data['substrate_strain'], start_coordinates=x.tolist(),
            start_norm_A=float(np.linalg.norm(x)), constraint_offset_A=expected.tolist(),
            held_values_A=held.tolist(), registration=registration,
            numerical_hessian_39=prep['numerical_hessian_39'])
        save(folder/'condition.json', row); write(folder/'start.traj', start)
        rows.append({k:v for k,v in row.items() if k not in ('numerical_hessian_39','start_coordinates','registration')})
    result = dict(policy=POLICY, diagonals=diagonal, conditions=rows,
        source_files={str(p.relative_to(REPO)).replace('\\','/'):sha256(p) for p in
            (OLD/'native_terminal.tar.gz',NEW/'completed/minus/native_terminal.tar.gz')},
        new_DFT_calls=0, candidate_energy_not_assigned_to_transported_starts=True)
    save(output/'config.json', result)
    return result


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); r=prepare(a.output); print(json.dumps(r['conditions']))
