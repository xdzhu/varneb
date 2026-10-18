"""Replay a completed F/A/J arm from portable native output, without DFT.

An unavailable arm has no candidate energy. A successful arm is a local
gradient candidate, not a revised path barrier or a certified minimum.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import sys

REPOSITORY = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPOSITORY))
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from examples.hfo2_fixed_input_factory import CONTRACT, read_fixed_hfo2_stru, same_ordered_geometry
from scripts.audit_hfo2_static_replica import audited_results, sha256


def review(case, exported, index, arm):
    case, exported = Path(case), Path(exported)
    if index not in (0, 1, 2) or arm not in ('A', 'J'):
        raise ValueError('one original T point and one declared arm required')
    manifest_path = case/'executed_bundle_files.json'
    if sha256(manifest_path) != json.loads((case/'package_receipt.json').read_bytes())['bundle_pin']:
        raise ValueError('executed bundle changed')
    files = json.loads(manifest_path.read_bytes())
    seed = case/'prepared/seeds'/str(index)
    if sha256(seed/'seed.json') != files[f'seeds/{index}/seed.json'] or sha256(case/'study.py') != files['study.py']:
        raise ValueError('executed input or source changed')
    spec = importlib.util.spec_from_file_location('E270_arm_data', case/'study.py')
    study = importlib.util.module_from_spec(spec); spec.loader.exec_module(study)
    metadata = json.loads((seed/'seed.json').read_bytes())
    path = exported/f'treatment_{arm}.json'; row = json.loads(path.read_bytes())
    if row['treatment'] != arm or not 0 <= row['new_attempts'] <= 12:
        raise ValueError('original treatment identity or budget differs')
    result = dict(index=index, strain=metadata['substrate_strain'], arm=arm,
        status=row['status'], total_arm_attempts_reported=row['new_attempts'],
        treatment_sha256=sha256(path), energy_change_meV_fu=None,
        released_gradient_norm_eV_A=None, local_native_SCFs_reparsed=0,
        new_DFT_calls=0, changed_path_label_or_forecast_score=False,
        energy_partition_or_stationary_TS_or_minimum_certified=False)
    if row['candidate'] is None:
        if row['status'] not in ('attempt_budget_exhausted', 'released_gradient_target_not_reached',
                                 'not_attempted_atomic_candidate_unavailable'):
            raise ValueError('unknown unavailable-candidate reason')
        return result
    if row['status'] != 'gradient_candidate_stability_unchecked':
        raise ValueError('only an original gradient candidate may have an energy')
    candidate = row['candidate']; x = np.asarray(candidate['coordinates'])
    chart = study.chart_for(seed); e0, g0 = study.cached_frozen(seed, chart)
    if x.shape != (39,) or not np.isfinite(x).all() or arm == 'A' and np.any(x[36:]):
        raise ValueError('finite fixed-chart coordinate and exact A cell required')
    if not np.any(x):
        energy, gradient = e0, g0
    else:
        matches = [p for p in exported.glob('evaluations/eval-*/result.json')
                   if json.loads(p.read_bytes())['coordinates'] == candidate['coordinates']]
        if len(matches) != 1:
            raise ValueError('candidate must match one complete native record')
        native_path = matches[0]; native = json.loads(native_path.read_bytes())
        raw = native_path.parent/'scf_000000'; validation = native['validation']
        if (native['status'] != 'complete' or validation['genuine_MPI_ranks'] != 32
                or any(validation['input_sha256'][name] != h for name, h in CONTRACT.items())
                or any(sha256(raw/name) != validation['input_sha256'][name] for name in ('INPUT', 'KPT', 'STRU'))
                or sha256(raw/'OUT.ABACUS/running_scf.log') != validation['raw_log_sha256']
                or sha256(raw/'call_audit.json') != validation['call_audit_sha256']):
            raise ValueError('portable native output/input identity changed')
        atoms = chart.to_atoms(x)
        if not same_ordered_geometry(atoms, read_fixed_hfo2_stru(raw/'STRU')):
            raise ValueError('accepted coordinate differs from its actual ordered input')
        atoms.calc = SinglePointCalculator(atoms, **audited_results(raw))
        energy, gradient = chart.energy_and_gradient(x, atoms)
        if (energy != native['enthalpy_eV']
                or not np.allclose(gradient, native['gradient'], atol=2e-12, rtol=0)):
            raise ValueError('native candidate E/physical-gradient replay differs')
        result.update(local_native_SCFs_reparsed=1, native_result_sha256=sha256(native_path),
                      final_native_call_seconds=validation['elapsed_seconds'],
                      full_orbital_pseudopotential_bytes_rechecked_locally=False)
    c = np.asarray(metadata['covectors_39_by_4']); q = np.asarray(metadata['orthonormal_39_by_4'])
    translation = np.zeros((39, 3)); translation[:36] = np.tile(np.eye(3), (12, 1))/np.sqrt(12)
    n = 36 if arm == 'A' else 39
    g = gradient[:n]
    released = g-q[:n]@(q[:n].T@g)-translation[:n]@(translation[:n].T@g)
    residual = float(np.linalg.norm(released)); defect = c.T@x
    if (energy != candidate['energy_eV_cell'] or residual > .03
            or not np.isclose(residual, candidate['released_gradient_norm_eV_A'], atol=2e-12, rtol=0)
            or np.max(np.abs(defect)) > 1e-10
            or candidate['stationary_TS_or_stability_certified'] is not False):
        raise ValueError('candidate energy, release target, fixed variables or claim differs')
    result.update(frozen_energy_eV_cell=e0, candidate_energy_eV_cell=energy,
        energy_change_meV_fu=250*(energy-e0), released_gradient_norm_eV_A=residual,
        physical_retained_defect_A=float(np.max(np.abs(defect))),
        coordinate_norm_A=float(np.linalg.norm(x)), cell_coordinate_norm_A=float(np.linalg.norm(x[36:])),
        complete_allowed_gradient_norm_eV_A=float(np.linalg.norm(gradient)),
        total_arm_cost_not_reconstructed_here=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exported', type=Path, required=True)
    parser.add_argument('--index', type=int, required=True, choices=(0, 1, 2))
    parser.add_argument('--arm', required=True, choices=('A', 'J'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = review(Path(__file__).resolve().parent, args.exported, args.index, args.arm)
    with args.output.open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(result, handle, indent=2, allow_nan=False); handle.write('\n')
    print(json.dumps(result))
