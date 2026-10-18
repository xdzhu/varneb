"""Read-only complete E272 native export, including budget-failed attempts."""
from __future__ import annotations
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import tarfile

import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from scripts.audit_hfo2_static_replica import audited_results, sha256

ROOT = Path('/public/home/iai806/abacus/agent-runs/20261018-varneb-preserving-branches-E272-r1')
PIN = 'd12af56eac84edf28e536c69b524d12681f876db1918c889fbf4d89eb59e06c2'
JOB = '28948425'


def save(path, value):
    with path.open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(value, handle, indent=2, allow_nan=False); handle.write('\n')


def terminal_row(raw):
    rows = [r.split('|') for r in raw.strip().splitlines() if r.split('|')[0] == JOB]
    if (len(rows) != 1 or len(rows[0]) < 8
            or rows[0][1:7] != ['COMPLETED', '0:0', '32', 'hfacnormal01', 'iai806', '1']):
        raise ValueError('only exact completed original E272 one-node32CPU job admitted')
    return rows[0]


def collect():
    raw = subprocess.check_output(['sacct', '-n', '-X', '-P', '-j', JOB,
        '--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,User,AllocNodes,ElapsedRaw,Start,End'], text=True)
    accounting = terminal_row(raw)
    spec = importlib.util.spec_from_file_location('E272_original_worker_export', ROOT/'bundle/worker.py')
    worker = importlib.util.module_from_spec(spec); spec.loader.exec_module(worker)
    study, _, config = worker.validate(PIN)
    complete = json.loads((ROOT/'results/complete.json').read_bytes())
    if complete['job_id'] != JOB or complete['bundle_pin'] != PIN:
        raise ValueError('original whole-batch terminal marker required')
    output = ROOT/'native_terminal'; output.mkdir(exist_ok=False)
    branches, costs = [], []
    for index in (0, 2):
        seed = ROOT/'bundle/prepared/seeds'/str(index)
        chart = study.chart_for(seed); e0, _ = study.cached_frozen(seed, chart)
        meta = json.loads((seed/'seed.json').read_bytes())
        q = np.asarray(meta['orthonormal_39_by_4']); c = np.asarray(meta['covectors_39_by_4'])
        translation = np.zeros((39, 3)); translation[:36] = np.tile(np.eye(3), (12, 1))/np.sqrt(12)
        for sign in ('minus', 'plus'):
            folder = ROOT/'results'/str(index)/sign; target = output/str(index)/sign
            target.mkdir(parents=True)
            summary = json.loads((folder/'treatments.json').read_bytes()); records = []
            for path in sorted(folder.glob('evaluations/eval-*/result.json')):
                native = json.loads(path.read_bytes()); atoms = chart.to_atoms(native['coordinates'])
                atoms.calc = SinglePointCalculator(atoms, **audited_results(path.parent/'scf_000000'))
                validation = study.result_validator(0, atoms, path.parent)
                energy, gradient = chart.energy_and_gradient(native['coordinates'], atoms)
                if (native['status'] != 'complete' or validation != native['validation']
                        or energy != native['enthalpy_eV']
                        or not np.allclose(gradient, native['gradient'], atol=2e-12, rtol=0)):
                    raise ValueError('native E/F/stress, six assets or physical-gradient replay differs')
                x = np.asarray(native['coordinates']); released = gradient[:36]-q[:36]@(q[:36].T@gradient[:36])
                released -= translation[:36]@(translation[:36].T@gradient[:36])
                record = dict(coordinates=native['coordinates'], energy_eV_cell=energy,
                    gradient=gradient.tolist(), native_result_sha256=sha256(path),
                    relative_record=path.relative_to(ROOT/'results').as_posix(),
                    coordinate_norm_A=float(np.linalg.norm(x)), atomic_released_gradient_eV_A=float(np.linalg.norm(released)))
                records.append(record); costs.append(validation['elapsed_seconds'])
                destination = target/path.parent.relative_to(folder)
                (destination/'scf_000000/OUT.ABACUS').mkdir(parents=True)
                for name in ('result.json', 'scf_000000/INPUT', 'scf_000000/KPT', 'scf_000000/STRU',
                             'scf_000000/call_audit.json', 'scf_000000/OUT.ABACUS/running_scf.log',
                             'scf_000000/abacus.out', 'scf_000000/abacus.err'):
                    if (path.parent/name).is_file(): shutil.copyfile(path.parent/name, destination/name)
            check = summary['known_direction_check']
            if (summary['index'] != index or summary['sign'] != sign or summary['strain'] != meta['substrate_strain']
                    or [r['treatment'] for r in summary['treatments']] != ['A', 'J']
                    or len(records) != summary['total_new_attempts']
                    or len(records) != sum(r['new_attempts'] for r in summary['treatments'])+check['new_attempts']
                    or len(records) > 26):
                raise ValueError('all branch attempts and finite native counts required')
            arms = []
            for arm, row in zip(('A', 'J'), summary['treatments']):
                candidate = row['candidate']; energy = residual = None
                if candidate is None and row['status'] not in ('attempt_budget_exhausted',
                        'released_gradient_target_not_reached', 'not_attempted_atomic_candidate_unavailable'):
                    raise ValueError('unknown unavailable-candidate reason')
                if candidate is not None:
                    matches = [r for r in records if r['coordinates'] == candidate['coordinates']]
                    if len(matches) != 1: raise ValueError('one exact native candidate required')
                    record = matches[0]; g = np.asarray(record['gradient']); x = np.asarray(record['coordinates'])
                    n = 36 if arm == 'A' else 39
                    residual = float(np.linalg.norm(g[:n]-q[:n]@(q[:n].T@g[:n])-translation[:n]@(translation[:n].T@g[:n])))
                    energy = record['energy_eV_cell']
                    if (row['status'] != 'gradient_candidate_stability_unchecked' or residual > .03
                            or energy != candidate['energy_eV_cell'] or np.max(np.abs(c.T@x)) > 1e-10
                            or not np.isclose(residual, candidate['released_gradient_norm_eV_A'], atol=2e-12, rtol=0)
                            or arm == 'A' and np.any(x[36:])):
                        raise ValueError('candidate target, held observables or exact A cell differs')
                arms.append(dict(arm=arm, status=row['status'], new_attempts=row['new_attempts'],
                    energy_eV_cell=energy, dE_minus_F_meV_fu=None if energy is None else 250*(energy-e0),
                    released_gradient_eV_A=residual))
            for name in ('treatments.json','treatment_A.json','treatment_J.json','known_direction_check.json'):
                shutil.copyfile(folder/name,target/name)
            branches.append(dict(index=index, strain=meta['substrate_strain'], sign=sign, arms=arms,
                known_direction_check=check, attempts=records, frozen_energy_eV_cell=e0))
    if len(costs) != complete['total_new_attempts'] or len(costs) > 104:
        raise ValueError('whole batch accounting differs')
    audit = dict(job_id=JOB,bundle_pin=PIN,accounting=accounting,branches=branches,
        new_native_SCF_calls=len(costs),native_call_proxy_core_hours=sum(costs)*32/3600,
        allocation_core_hours=float(accounting[7])*32/3600, all_attempts_exported=True,
        all_native_six_assets_MPI_EFS_and_chart_gradients_replayed_on_HF=True,
        old_preparation_and_curvature_costs_additional=True, new_DFT_calls=0,scheduler_mutations=0,
        original_labels_and_forecast_scores_unchanged=True, stationary_TS_or_full_stability_certified=False)
    save(output/'audit.json',audit); shutil.copyfile(ROOT/'results/complete.json',output/'complete.json')
    save(output/'files.json',{p.relative_to(output).as_posix():sha256(p) for p in sorted(output.rglob('*')) if p.is_file()})
    archive=ROOT/'native_terminal.tar.gz'
    with tarfile.open(archive,'x:gz') as tar: tar.add(output,arcname='native_terminal')
    receipt=dict(archive_sha256=sha256(archive), archive_bytes=archive.stat().st_size,
        exported_files=len(json.loads((output/'files.json').read_bytes())),new_SCF_calls=len(costs),
        allocation_core_hours=audit['allocation_core_hours'],native_call_proxy_core_hours=audit['native_call_proxy_core_hours'],
        job_id=JOB,new_DFT_calls=0,all_attempts_exported=True)
    save(ROOT/'native_terminal_receipt.json',receipt)
    return receipt


if __name__ == '__main__': print(json.dumps(collect()))
