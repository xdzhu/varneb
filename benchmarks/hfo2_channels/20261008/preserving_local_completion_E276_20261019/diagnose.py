"""Replay the two original +1% failed arms; no DFT, fit or resubmission."""
from __future__ import annotations
import argparse
import csv
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import numpy as np

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))
from scripts.audit_hfo2_static_replica import sha256

CASE = Path(__file__).resolve().parent
OLD = CASE.with_name('preserving_pattern_branches_E272_20261018')


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


def analyze():
    reader = module(OLD/'analyze_prefix.py', 'E276_native_reader')
    terminal = module(OLD/'analyze_terminal.py', 'E276_archive_reader')
    chart, data, e0, _ = reader.setup(2)
    q = np.asarray(data['orthonormal_39_by_4'])
    t = np.zeros((39, 3)); t[:36] = np.tile(np.eye(3), (12, 1))/np.sqrt(12)
    project = np.eye(39)-q@q.T-t@t.T; project[36:] = 0.; project[:, 36:] = 0.
    covectors = np.asarray(data['covectors_39_by_4'])
    rows, arms = [], []
    with tempfile.TemporaryDirectory(prefix='E276-original-failures-') as td:
        root, receipt = terminal.restore(Path(td))
        for sign in ('minus', 'plus'):
            folder = root/'2'/sign
            files = sorted(folder.glob('evaluations/eval-*/result.json'))
            if len(files) != 12: raise ValueError('the exact two original twelve-call arms are required')
            previous = None; trace = []
            for path in files:
                native = json.loads(path.read_bytes()); x = np.asarray(native['coordinates'])
                atoms, energy, gradient, origin = reader.replay_native(folder, chart, x.tolist())
                if np.any(x[36:]) or np.max(np.abs(covectors.T@x)) > 1e-10:
                    raise ValueError('original atomic-only fixed observables changed')
                g = project@gradient; norm = float(np.linalg.norm(x)); radial = x/norm
                radial_g = float(g@radial)
                row = dict(sign=sign, evaluation=path.parent.name, norm_A=norm,
                    dE_minus_F_meV_fu=250*(energy-e0), released_gradient_eV_A=float(np.linalg.norm(g)),
                    radial_gradient_eV_A=radial_g,
                    outward_downhill_component_eV_A=-radial_g,
                    tangent_to_sphere_gradient_eV_A=float(np.linalg.norm(g-radial_g*radial)),
                    step_from_previous_trial_A=None if previous is None else float(np.linalg.norm(x-previous)),
                    retained_defect_A=float(np.max(np.abs(covectors.T@x))),
                    source_result_sha256=origin['native_result_sha256'])
                rows.append(row); trace.append(row); previous = x
            last = json.loads(files[-1].read_bytes())
            arms.append(dict(sign=sign, last_coordinates=last['coordinates'],
                last_energy_eV_cell=last['enthalpy_eV'], last_gradient=last['gradient'],
                trace=trace, last_native_relative_path=str(files[-1].relative_to(root)).replace('\\','/'),
                last_result_sha256=sha256(files[-1]), final_is_not_a_qualified_candidate=True))
    return dict(status='two_original_failed_atomic_arms_replayed_without_DFT',
        original_archive_sha256=receipt['archive_sha256'], original_source_job='28812952',
        original_local_job='28948425', original_norm_guard_A=.5,
        original_fresh_attempt_cap=12, original_gradient_tolerance_eV_A=.03,
        original_physical_input_sha256=data['physical_input_sha256'],
        new_DFT_calls=0, geometry_or_physical_constraints_changed=False,
        traces_are_evaluated_trials_not_accepted_optimizer_iterates=True,
        sphere_tangent_gradient_not_used_as_original_convergence_criterion=True,
        original_F_energy_eV_cell=e0, local_native_calls_reparsed=24, arms=arms), rows


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    args=p.parse_args(); result, rows=analyze(); args.output.mkdir(parents=True, exist_ok=False)
    with (args.output/'diagnosis.json').open('x',encoding='utf-8',newline='\n') as f:
        json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    with (args.output/'trial_trace.csv').open('x',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    for arm in result['arms']:
        print(json.dumps(dict(sign=arm['sign'], first=arm['trace'][0], last=arm['trace'][-1])))
