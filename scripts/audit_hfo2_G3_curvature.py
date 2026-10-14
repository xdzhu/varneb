"""Read-only raw-output audit of the fixed E150 local-curvature stencil.

A live complete-result prefix may be audited, but missing columns are never
filled. No calculator, scheduler operation, reference selection or TS claim.
Run against the unchanged E150 runtime through PYTHONPATH; write only a fresh
external report, never inside the running study or its source directories.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator

import scripts.hfo2_G3_local_curvature as study
from vcneb import assemble_joint_directional_curvature
from vcneb.mode_surface import (audit_directional_curvature_consistency,
                                audit_directional_work_consistency)


def _plain(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(type(value).__name__)


def diagnostics(manifest, channel, center_energy, center_gradient, rows):
    """Numeric diagnostics of an already raw-audited, ordered complete prefix.

    This function alone does not establish provenance. Per-step incomplete
    projected blocks are deliberately not assigned a full internal index.
    """
    _, probes = study.validate(manifest, channel)
    n = len(rows)
    gradient = np.asarray(center_gradient, dtype=float)
    if (n > 148 or gradient.shape != (40,) or not np.isfinite(gradient).all()
            or not np.isfinite(center_energy)):
        raise ValueError('finite original center and at most148 records required')
    for index, row in enumerate(rows):
        values = np.asarray(row['coordinates'], dtype=float)
        g = np.asarray(row['gradient'], dtype=float)
        if (row['status'] != 'complete' or values.shape != (40,)
                or not np.array_equal(values, probes[index].delta_A)
                or g.shape != (40,) or not np.isfinite(g).all()
                or not np.isfinite(row['enthalpy_eV'])):
            raise ValueError('ordered complete prefix differs from registered stencil')
    basis = np.asarray(manifest['centers'][channel]['directions'])
    blocks = []
    for step_index, h in enumerate(study.STEPS):
        start = step_index * 74
        pairs = max(0, min(74, n-start)) // 2
        if not pairs:
            continue
        g = np.asarray([r['gradient'] for r in rows[start:start+2*pairs]])
        block = assemble_joint_directional_curvature(basis[:, :pairs], g[::2], g[1::2], step_A=h)
        energies = np.asarray([r['enthalpy_eV'] for r in rows[start:start+2*pairs]]).reshape(pairs, 2)
        energy_diagonal = (energies.sum(axis=1)-2*center_energy)/h**2
        complete = pairs == 37
        blocks.append(dict(step_A=h, complete_pairs=pairs, full_stencil_at_this_step=complete,
            raw_projected_eV_A2=block.raw_projected,
            full_hessian_action_eV_A2=block.full_hessian_action,
            reciprocity_relative_defect=block.reciprocity_relative_defect,
            symmetrization_operator_change_eV_A2=block.symmetrization_operator_change_eV_A2,
            energy_diagonal_eV_A2=energy_diagonal,
            energy_gradient_diagonal_defect_eV_A2=energy_diagonal-np.diag(block.raw_projected),
            fixed_strain_internal_eigenvalues_eV_A2=(
                np.linalg.eigvalsh(block.symmetric_projected[:36, :36]) if complete else None)))
    lines = []
    for j in range(max(0, min(74, n-74)) // 2):
        indices = [2*j, 2*j+1, 74+2*j, 75+2*j]
        e = np.asarray([rows[i]['enthalpy_eV'] for i in indices]).reshape(2, 2)
        g = np.asarray([np.asarray(rows[i]['gradient']) @ basis[:, j] for i in indices]).reshape(2, 2)
        consistency = audit_directional_curvature_consistency(center_energy, e, g, study.STEPS)
        work = audit_directional_work_consistency(
            np.array([-.02, -.01, 0., .01, .02]),
            np.array([e[1, 1], e[0, 1], center_energy, e[0, 0], e[1, 0]]),
            np.array([g[1, 1], g[0, 1], gradient @ basis[:, j], g[0, 0], g[1, 0]]))
        lines.append(dict(direction=j, control_direction=j == 36,
            energy_gradient_curvature=asdict(consistency), five_point_work=asdict(work)))
    comparison = None
    if n == 148:
        matrices = [np.asarray(b['raw_projected_eV_A2']) for b in blocks]
        internal = [(m[:36, :36]+m[:36, :36].T)/2 for m in matrices]
        tangent = basis[:, :36].T @ np.asarray(manifest['centers'][channel]['ordered_local_tangent'])
        values, vectors = zip(*(np.linalg.eigh(m) for m in internal))
        comparison = dict(
            internal_operator_step_change_eV_A2=float(np.linalg.norm(internal[1]-internal[0], ord=2)),
            raw_full_projected_step_change_eV_A2=float(np.linalg.norm(matrices[1]-matrices[0], ord=2)),
            lowest_internal_direction_absolute_overlap=float(abs(vectors[0][:, 0] @ vectors[1][:, 0])),
            lowest_direction_native_tangent_absolute_overlap=[float(abs(v[:, 0] @ tangent)) for v in vectors],
            internal_counts_below_zero=[int(np.count_nonzero(v < 0)) for v in values],
            counts_are_raw_not_resolution_certified=True, external_control_excluded=True)
    result = dict(status='complete_local_curvature_diagnostic' if n == 148 else 'partial_local_curvature_diagnostic',
        completed_records=n, planned_records=148, blocks=blocks, two_step_lines=lines,
        two_step_comparison=comparison, raw_center_gradient_eV_A=gradient,
        physical_error_bound_or_TS_certificate=False, selected_conditional_plane=False,
        limitations=['Incomplete projected blocks do not count full-space instabilities.',
            'Raw reciprocity, step spread and energy/force mismatch are diagnostics, not total DFT error bounds.',
            'Nonstationary local curvature and tangent overlap do not certify a TS or basin connections.',
            'No unmeasured translation self-columns or ASR repair; no mode/frequency identification.'])
    return json.loads(json.dumps(result, default=_plain, allow_nan=False))


def audit(manifest, channel, workdir):
    """Reparse one immutable native result prefix; never instantiate a calculator."""
    chart, probes, energy, gradient = study.native_preflight(manifest, channel)
    center_path = workdir/'center.json'
    center = json.loads(center_path.read_text())
    if (center['cached_energy_eV_cell'] != energy or center['new_center_SCFs'] != 0
            or not np.array_equal(center['cached_full_gradient_eV_A'], gradient)):
        raise ValueError('cached native center differs')
    directory = workdir/'evaluations'
    contract = json.loads((directory/'contract.json').read_text())
    fingerprint = hashlib.sha256(json.dumps({k:v for k,v in contract.items() if k != 'contract_sha256'},
        sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    if (contract.get('contract_sha256') != fingerprint
            or contract.get('coordinate_convention') != chart.coordinate_convention
            or contract.get('mechanical_boundary') != chart.mechanical_contract
            or contract.get('pressure_eV_per_A3') != 0.
            or contract.get('calculator_id') != 'original-six-HfO2-PBE-100Ry-full10auDZP-Gamma2x2x2'
            or contract.get('validation_id') != 'E150-native-six-input-DSIZE32-EFS-one-SCF-v1'):
        raise ValueError('native evaluator contract differs')
    paths = sorted(directory.glob('eval-*/result.json'))  # atomic complete files only
    rows, sources = [], []
    for index, path in enumerate(paths):
        if index >= 148:
            raise ValueError('more than148 completed records; do not silently truncate')
        row = json.loads(path.read_text())
        digest = hashlib.sha256(probes[index].delta_A.tobytes()).hexdigest()
        if (path.parent.name != f'eval-{index:06d}-{digest[:12]}'
                or row['coordinate_sha256'] != digest
                or row['contract_sha256'] != fingerprint):
            raise ValueError('missing, duplicated or remapped native prefix record')
        raw = study.audited_results(path.parent/'scf_000000')
        atoms = probes[index].atoms
        atoms.calc = SinglePointCalculator(atoms, **raw)
        validation = study.result_validator(index, atoms, path.parent)
        e, g = chart.energy_and_gradient(probes[index].delta_A, atoms)
        if (validation != row['validation'] or e != row['enthalpy_eV']
                or not np.array_equal(g, row['gradient'])
                or not np.array_equal(atoms.get_stress(voigt=False), row['stress_eV_per_A3'])):
            raise ValueError('cached result differs from reparsed native E/F/stress/gradient')
        rows.append(row)
        sources.append(dict(result_sha256=study.sha256(path), validation=validation))
    result = diagnostics(manifest, channel, energy, gradient, rows)
    failures = [str(p) for p in directory.glob('eval-*/failure.json')]
    if failures:
        result['status'] = 'native_evaluation_failure_present'
    result.update(channel=channel, batch_job_id=center['job_id'], raw_SCFs_reparsed=len(rows),
        native_record_sources=sources, center_sha256=study.sha256(center_path),
        center_source=manifest['centers'][channel]['raw_source'],
        evaluator_contract_sha256=fingerprint, failures_present=failures,
        conservative_attempt_reservations=len(list(directory.glob('attempt-*.json'))),
        observed_native_SCF_directories=len(list(directory.glob('eval-*/scf_*'))),
        sum_audited_call_seconds=sum(r['validation']['elapsed_seconds'] for r in rows),
        new_DFT_or_scheduler_operations=0, allocation_core_hours=None,
        physical_reserved_geometry_generated_or_read=False)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--manifest-sha256', required=True)
    parser.add_argument('--channel', choices=tuple(study.SOURCES), required=True)
    parser.add_argument('--workdir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('fresh audit report required')
    if any(args.output.resolve().is_relative_to(p.resolve()) for p in
            (args.workdir.parent, study.ASSETS)):
        raise ValueError('never write an audit inside the live study/source')
    if study.sha256(args.manifest) != args.manifest_sha256:
        raise ValueError('frozen manifest byte identity differs')
    result = audit(json.loads(args.manifest.read_text()), args.channel, args.workdir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    study.exclusive_json(args.output, result)
    print(json.dumps({k:result[k] for k in ('status', 'channel', 'raw_SCFs_reparsed', 'failures_present')}))


if __name__ == '__main__':
    main()
