"""Terminal E268 native audit and sampled-value scoring, never a TS certificate."""
import argparse
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess

import numpy as np

from examples.hfo2_fixed_input_factory import CONTRACT
from scripts.audit_hfo2_static_replica import audited_results, sha256
from scripts.audit_hfo2_clamped_terminal import completed_step, validate_terminal
from scripts.export_hfo2_clamped_observation import export, REGISTERED_FACTORIES
from vcneb.forecast_freeze import load_frozen_forecast_panel

ROOT = Path('/public/home/iai806/abacus/agent-runs/20261018-varneb-physical-E268-r1')
CHANNELS = ('PO_flip_T_pattern_preserving', 'PO_to_T')


def save(path, value):
    with path.open('x', encoding='utf-8') as f:
        json.dump(value, f, indent=2, allow_nan=False); f.write('\n')


def audit(index, output):
    if output.exists(): raise FileExistsError(output)
    accepted = json.loads((ROOT/f'path_accepted_{index}.json').read_bytes())
    pin = accepted['path_bundle_sha256']
    control_file = ROOT/'path_bundle/path_control.py'
    files_path = ROOT/'path_bundle/bundle_files.json'
    if sha256(files_path) != pin or sha256(control_file) != json.loads(files_path.read_bytes())['path_control.py']:
        raise ValueError('immutable production control changed')
    spec = importlib.util.spec_from_file_location('E268_frozen_control', control_file)
    control = importlib.util.module_from_spec(spec); spec.loader.exec_module(control)
    config, _ = control.validate(pin)
    job = accepted['job_id']; channel = CHANNELS[index]
    if accepted['channel'] != channel: raise ValueError('original channel handle mismatch')
    raw = subprocess.check_output(['sacct', '-n', '-X', '-P', '-j', job,
        '--format=JobIDRaw,State,ExitCode,AllocCPUS,Partition,User,AllocNodes,ElapsedRaw,Start,End'], text=True).strip()
    row = raw.split('|')
    if len(row) < 10 or row[6] != '1': raise ValueError('complete one-node accounting row required')
    work = ROOT/'paths'/channel/'band'
    summary = json.loads((work/'vcneb_summary.json').read_bytes())
    outcome = 'converged' if summary.get('converged') else 'step-cap'
    validate_terminal(job, '|'.join(row[:6]), summary, (work/'vcneb_failure.json').exists(), outcome=outcome)
    step = completed_step(work/'vcneb.opt.log')
    if step > 80: raise ValueError('registered80-update cap exceeded')
    calls = sorted(work.glob('image_*/scf_*/call_audit.json'))
    if len(calls) != 7*(step+1): raise ValueError('fresh interior evaluation count differs')
    call_records = []
    for path in calls:
        a = json.loads(path.read_bytes()); folder = path.parent
        if (any(a['input_sha256'][n] != h or sha256(folder/n) != h for n, h in CONTRACT.items())
                or sha256(folder/'STRU') != a['input_sha256']['STRU']
                or sha256(folder/'OUT.ABACUS/running_scf.log') != a['raw_log_sha256']):
            raise ValueError('native input/log hash drift')
        results = audited_results(folder)
        if any(not np.allclose(results[k], a['results'][k], atol=1e-12, rtol=0) for k in results):
            raise ValueError('native E/F/stress differs')
        if not np.isfinite(a['elapsed_seconds']) or a['elapsed_seconds'] <= 0: raise ValueError('invalid call cost')
        call_records.append(dict(directory=str(folder), call_audit_sha256=sha256(path),
            input_sha256=a['input_sha256'], raw_log_sha256=a['raw_log_sha256'], elapsed_seconds=a['elapsed_seconds']))
    script = ROOT/'path_bundle/paths.slurm'
    # In-memory registration of this separately frozen new recipe only.
    # No old script, old accepted hash or source file is changed.
    script_pin = sha256(script)
    REGISTERED_FACTORIES[script_pin] = 'examples.hfo2_fixed_input_factory:make_clamped_seed_cached_factory'
    output.mkdir()
    observation = export(work, step, output/f'step_{step:04d}', job, production_script=script)
    if bool(observation['ordinary_residual_pass']) != bool(summary['converged']):
        raise ValueError('ordinary physical replay and summary differ')
    for item in observation['raw_image_evaluations']:
        source = Path(item['raw_source']); dest = output/f'step_{step:04d}'/'raw'/f"image_{item['image_index']:04d}"
        (dest/'OUT.ABACUS').mkdir(parents=True)
        for name in ('INPUT', 'KPT', 'STRU', 'OUT.ABACUS/running_scf.log'):
            shutil.copyfile(source/name, dest/name)
    shutil.copyfile(work/'vcneb_summary.json', output/'vcneb_summary.json')
    result = dict(status='audited_ordinary_sampled_holdout' if summary['converged'] else 'audited_bounded_step_cap',
        job_id=job, channel=channel, scheduler_row=row, step=step, ordinary_converged=summary['converged'],
        final_fmax_eV_A=observation['replayed_fmax_eV_A'],
        sampled_peak_eV_fu=observation['sampled_forward_barrier_meV_fu']/1000,
        physical_error_bound_eV_fu=None, TS_certified=False, new_DFT_calls=0,
        native_all_calls_reparsed=True, all_six_assets_checked=True, real_MPI32_checked=True,
        calls=call_records, new_SCF_calls=len(calls),
        SCF_core_hours=sum(r['elapsed_seconds'] for r in call_records)*32/3600,
        allocation_core_hours=float(row[7])*32/3600, prior_endpoint_and_model_costs_additional=True,
        path_bundle_sha256=pin, forecast_manifest_sha256=config['forecast_manifest_sha256'],
        observation_sha256=sha256(output/f'step_{step:04d}'/'observation.json'))
    save(output/'audit.json', result); return result


def point_scores(forecasts, references):
    """Descriptive point errors only; absent uncertainty stays absent."""
    expected = set(CHANNELS)
    if set(references) != expected: raise ValueError('both registered outcomes required')
    models = sorted({r.model_id for r in forecasts}); rows = []
    for model in models:
        model_rows = [r for r in forecasts if r.model_id == model]
        if len(model_rows) != 2 or {r.case_id for r in model_rows} != expected:
            raise ValueError('complete forecast coverage including abstentions required')
        errors, points = [], {}
        for r in model_rows:
            ref = references[r.case_id]
            valid = ref['ordinary_converged'] and r.prediction_eV_fu is not None and r.abstention_reason is None
            error = (r.prediction_eV_fu-ref['sampled_peak_eV_fu'])*1000 if valid else None
            if valid: errors.append(abs(error))
            points[r.case_id] = dict(prediction_eV_fu=r.prediction_eV_fu,
                reference_sampled_peak_eV_fu=ref['sampled_peak_eV_fu'],
                reference_ordinary_converged=ref['ordinary_converged'],
                abstention_reason=r.abstention_reason, signed_error_meV_fu=error,
                response_error_meV_fu=error, physical_error_bound_eV_fu=None)
        gap_error = None
        if len(errors) == 2:
            gap_error = points['PO_to_T']['signed_error_meV_fu']-points['PO_flip_T_pattern_preserving']['signed_error_meV_fu']
        rows.append(dict(model=model, denominator=2, scored=len(errors), points=points,
            maximum_absolute_error_meV_fu=max(errors) if errors else None,
            mean_absolute_error_meV_fu=float(np.mean(errors)) if errors else None,
            T_minus_preserving_error_meV_fu=gap_error,
            physically_resolved_advantage_certified=False))
    return rows


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--index', type=int, choices=(0, 1), required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); r = audit(a.index, a.output)
    print(json.dumps({k: v for k, v in r.items() if k != 'calls'}))
