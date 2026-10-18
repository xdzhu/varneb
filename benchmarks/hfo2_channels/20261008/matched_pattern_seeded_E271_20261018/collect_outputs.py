"""Portable native export of one terminal E270/E271 batch, without DFT."""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import tarfile

import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from scripts.audit_hfo2_static_replica import audited_results, sha256

BASE = Path('/public/home/iai806/abacus/agent-runs/20261018-varneb-matched-pattern-E270-r1')
ROOT = BASE.with_name('20261018-varneb-matched-pattern-E271-r1')
PINS = dict(E270='f2aca3403e214fbad072ed53e69cd1afc496c560e52f2909aabe233dfb2b42ca',
            E271='bcdb328562237d84522ddc5aa597b291ef5ba622716907e6692cd530cbda32f6')


def save(path, value):
    with path.open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(value, handle, indent=2, allow_nan=False); handle.write('\n')


def require_terminal(raw, job):
    rows = [r.split('|') for r in raw.strip().splitlines() if r.split('|')[0] == job]
    if (len(rows) != 1 or len(rows[0]) < 8
            or rows[0][1:7] != ['COMPLETED', '0:0', '32', 'hfacnormal01', 'iai806', '1']):
        raise ValueError('exact terminal one-node32CPU job required; no live/timeout promotion')
    return rows[0]


def candidate_summary(study, chart, data, arm, row, candidates):
    if row['treatment'] != arm or not 0 <= row['new_attempts'] <= 12:
        raise ValueError('original treatment identity and finite attempt cap required')
    result = dict(arm=arm, status=row['status'], new_attempts=row['new_attempts'],
        candidate_energy_eV_cell=None, energy_change_from_F_meV_fu=None,
        released_gradient_eV_A=None, candidate_source=None)
    point = row['candidate']
    if point is None:
        if row['status'] not in ('attempt_budget_exhausted', 'released_gradient_target_not_reached',
                                 'not_attempted_atomic_candidate_unavailable'):
            raise ValueError('unknown unavailable-candidate status')
        return result
    if row['status'] != 'gradient_candidate_stability_unchecked':
        raise ValueError('successful local gradient candidate status required')
    x = np.asarray(point['coordinates']); matches = [r for r in candidates if r['coordinates'] == point['coordinates']]
    if x.shape != (39,) or not np.isfinite(x).all() or len(matches) != 1:
        raise ValueError('one exact complete native record must support each candidate')
    native = matches[0]; g = np.asarray(native['gradient'])
    q = np.asarray(data['orthonormal_39_by_4']); c = np.asarray(data['covectors_39_by_4'])
    t = np.zeros((39, 3)); t[:36] = np.tile(np.eye(3), (12, 1))/np.sqrt(12)
    n = 36 if arm == 'A' else 39
    residual = float(np.linalg.norm(g[:n]-q[:n]@(q[:n].T@g[:n])-t[:n]@(t[:n].T@g[:n])))
    if (point['energy_eV_cell'] != native['enthalpy_eV'] or residual > .03
            or not np.isclose(residual, point['released_gradient_norm_eV_A'], atol=2e-12, rtol=0)
            or np.max(np.abs(c.T@x)) > 1e-10 or np.max(np.abs(t.T@x)) > 1e-10
            or arm == 'A' and np.any(x[36:])
            or point['stationary_TS_or_stability_certified'] is not False):
        raise ValueError('local candidate energy, gradient or physical constraints differ')
    result.update(candidate_energy_eV_cell=native['enthalpy_eV'],
        energy_change_from_F_meV_fu=250*(native['enthalpy_eV']-data['frozen_energy_eV_cell']),
        released_gradient_eV_A=residual, physical_retained_defect_A=float(np.max(np.abs(c.T@x))),
        coordinate_norm_A=float(np.linalg.norm(x)), cell_coordinate_norm_A=float(np.linalg.norm(x[36:])),
        complete_allowed_gradient_eV_A=float(np.linalg.norm(g)), candidate_source=native['export_source'])
    return result


def collect(tag):
    root = BASE if tag == 'E270' else ROOT
    pin = PINS[tag]; bundle = root/'bundle'
    accepted = json.loads((root/'accepted.json').read_bytes()); job = accepted['job_id']
    if accepted['bundle_pin'] != pin or tag == 'E270' and job != '28947338':
        raise ValueError('only the accepted original handle admitted')
    raw = subprocess.check_output(['sacct', '-n', '-X', '-P', '-j', job,
        '--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,User,AllocNodes,ElapsedRaw,Start,End'], text=True)
    accounting = require_terminal(raw, job)
    complete = json.loads((root/'results/complete.json').read_bytes())
    if complete['job_id'] != job or complete['bundle_pin'] != pin:
        raise ValueError('terminal result/source identity differs')
    if sha256(bundle/'bundle_files.json') != pin:
        raise ValueError('executed manifest differs')
    for name, digest in json.loads((bundle/'bundle_files.json').read_bytes()).items():
        if Path(name).is_absolute() or '..' in Path(name).parts or sha256(bundle/name) != digest:
            raise ValueError('executed source/input bytes changed')
    spec = importlib.util.spec_from_file_location('immutable_E270_export', BASE/'bundle/study.py')
    study = importlib.util.module_from_spec(spec); spec.loader.exec_module(study)
    study.validate(PINS['E270'])
    output = root/'native_terminal'; output.mkdir(exist_ok=False)
    reports, calls = [], []
    for i in range(3):
        seed = BASE/'bundle/seeds'/str(i); chart = study.chart_for(seed)
        energy, frozen_gradient = study.cached_frozen(seed, chart)
        meta = json.loads((seed/'seed.json').read_bytes()); meta['frozen_energy_eV_cell'] = energy
        directory = root/'results'/str(i); destination = output/str(i); destination.mkdir()
        summary = json.loads((directory/'treatments.json').read_bytes())
        if summary['index'] != i or [r['treatment'] for r in summary['treatments']] != ['A', 'J']:
            raise ValueError('both ordered treatments at the original point required')
        records = []
        for path in sorted(directory.glob('evaluations/eval-*/result.json')):
            native = json.loads(path.read_bytes()); atoms = chart.to_atoms(native['coordinates'])
            atoms.calc = SinglePointCalculator(atoms, **audited_results(path.parent/'scf_000000'))
            validation = study.result_validator(0, atoms, path.parent)
            e, g = chart.energy_and_gradient(native['coordinates'], atoms)
            if (native['status'] != 'complete' or validation != native['validation']
                    or e != native['enthalpy_eV'] or not np.allclose(g, native['gradient'], atol=2e-12, rtol=0)):
                raise ValueError('native E/F/stress replay differs')
            relative = path.parent.relative_to(directory)
            target = destination/relative; (target/'scf_000000/OUT.ABACUS').mkdir(parents=True)
            for name in ('result.json', 'scf_000000/INPUT', 'scf_000000/KPT', 'scf_000000/STRU',
                         'scf_000000/call_audit.json', 'scf_000000/OUT.ABACUS/running_scf.log'):
                shutil.copyfile(path.parent/name, target/name)
            for name in ('abacus.out', 'abacus.err'):
                if (path.parent/'scf_000000'/name).is_file():
                    shutil.copyfile(path.parent/'scf_000000'/name, target/'scf_000000'/name)
            native['export_source'] = str(path); records.append(native)
            calls.append(dict(index=i, result_sha256=sha256(path), validation=validation))
        if (len(records) != summary['total_new_attempts']
                or len(records) != sum(r['new_attempts'] for r in summary['treatments'])):
            raise ValueError('all attempts, including numerical failures, must be exported')
        for name in ('treatments.json', 'treatment_A.json', 'treatment_J.json'):
            shutil.copyfile(directory/name, destination/name)
        arms = []
        for arm, row in zip(('A', 'J'), summary['treatments']):
            sources = records+[dict(coordinates=np.zeros(39).tolist(), enthalpy_eV=energy,
                gradient=frozen_gradient.tolist(), export_source=str(seed/'raw/STRU'))]
            if (tag == 'E271' and row['candidate'] is not None
                    and not any(r['coordinates'] == row['candidate']['coordinates'] for r in sources)):
                sources = []
                for p in (BASE/'results'/str(i)).glob('evaluations/eval-*/result.json'):
                    old = json.loads(p.read_bytes())
                    if old['coordinates'] == row['candidate']['coordinates']:
                        atoms = chart.to_atoms(old['coordinates'])
                        atoms.calc = SinglePointCalculator(atoms, **audited_results(p.parent/'scf_000000'))
                        if study.result_validator(0, atoms, p.parent) != old['validation']:
                            raise ValueError('reused original candidate changed')
                        e, g = chart.energy_and_gradient(old['coordinates'], atoms)
                        if e != old['enthalpy_eV'] or not np.allclose(g, old['gradient'], atol=2e-12, rtol=0):
                            raise ValueError('reused original gradient changed')
                        old['export_source'] = str(p); sources.append(old)
            arms.append(candidate_summary(study, chart, meta, arm, row, sources))
        reports.append(dict(index=i, strain=meta['substrate_strain'], frozen_energy_eV_cell=energy, arms=arms))
    shutil.copyfile(root/'results/complete.json', output/'complete.json')
    result = dict(tag=tag, job_id=job, bundle_pin=pin, accounting=accounting, rows=reports,
        new_native_SCF_calls=len(calls), native_call_proxy_core_hours=sum(r['validation']['elapsed_seconds'] for r in calls)*32/3600,
        allocation_core_hours=float(accounting[7])*32/3600, calls=calls,
        all_new_native_EFS_and_six_assets_checked_on_HF=True, old_costs_additional=True,
        old_path_labels_unchanged=True, stationary_TS_or_minimum_or_stable_surface_certified=False,
        new_DFT_calls=0, scheduler_mutations=0)
    save(output/'audit.json', result)
    files = {p.relative_to(output).as_posix(): sha256(p) for p in sorted(output.rglob('*')) if p.is_file()}
    save(output/'files.json', files)
    archive = root/'native_terminal.tar.gz'
    with tarfile.open(archive, 'x:gz') as tar: tar.add(output, arcname='native_terminal')
    receipt = dict(archive_sha256=sha256(archive), archive_bytes=archive.stat().st_size,
        exported_files=len(files), job_id=job, new_SCF_calls=len(calls),
        native_call_proxy_core_hours=result['native_call_proxy_core_hours'],
        allocation_core_hours=result['allocation_core_hours'], new_DFT_calls=0)
    save(root/'native_terminal_receipt.json', receipt)
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('tag', choices=('E270', 'E271'))
    print(json.dumps(collect(parser.parse_args().tag)))
