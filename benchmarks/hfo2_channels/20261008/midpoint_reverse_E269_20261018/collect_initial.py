"""Export the two complete step-zero snapshots without stopping live jobs.

These are explicitly unrelaxed reverse-strain seeds, not completed paths or
replacement holdout labels. The material runtime and its source stay unchanged.
"""
import importlib.util
import json
from pathlib import Path
import shutil
import tarfile

from collect_results import ROOT, PIN, JOBS, CHANNELS, accounting, save
from scripts.audit_hfo2_static_replica import sha256
from scripts.export_hfo2_clamped_observation import export, REGISTERED_FACTORIES


def collect():
    output, archive = ROOT/'initial_complete_chains', ROOT/'initial_complete_chains.tar.gz'
    if output.exists() or archive.exists():
        raise FileExistsError('initial export already exists; do not overwrite')
    spec = importlib.util.spec_from_file_location('E269_fixed_control', ROOT/'bundle/control.py')
    control = importlib.util.module_from_spec(spec); spec.loader.exec_module(control)
    control.validate(PIN)
    script = ROOT/'bundle/paths.slurm'
    REGISTERED_FACTORIES[sha256(script)] = 'examples.hfo2_fixed_input_factory:make_clamped_seed_cached_factory'
    rows = [accounting(job) for job in JOBS]
    for i, row in enumerate(rows):
        accepted = json.loads((ROOT/f'accepted_{i}.json').read_bytes())
        if (accepted['job_id'] != JOBS[i] or accepted['bundle_sha256'] != PIN
                or row[1] not in ('RUNNING', 'COMPLETED') or row[2] != '0:0'
                or row[3:7] != ['32', 'hfacnormal01', 'iai806', '1']):
            raise ValueError('original live or successful32CPU handle required')
        log = ROOT/'paths'/CHANNELS[i]/'band/vcneb.opt.log'
        if not any(line.split()[:2] == ['CheckedFIRE:', '0'] for line in log.read_text().splitlines()):
            raise ValueError('first complete evaluation is not yet recorded')
    output.mkdir()
    results = []
    for i, short in enumerate(('preserving', 'T')):
        work = ROOT/'paths'/CHANNELS[i]/'band'; folder = output/short; folder.mkdir()
        observation = export(work, 0, folder/'step_0000', JOBS[i], production_script=script)
        for image in observation['raw_image_evaluations']:
            source = Path(image['raw_source']); destination = folder/'step_0000/raw'/f"image_{image['image_index']:04d}"
            (destination/'OUT.ABACUS').mkdir(parents=True)
            for name in ('INPUT', 'KPT', 'STRU', 'OUT.ABACUS/running_scf.log'):
                shutil.copyfile(source/name, destination/name)
        report = dict(status='complete_initial_evaluation_not_terminal_result', channel=CHANNELS[i],
            job_id=JOBS[i], scheduler_row_at_export=rows[i], step=0, bundle_sha256=PIN,
            initializer_strain=.01, target_strain=.005, ordinary_residual_pass=observation['ordinary_residual_pass'],
            replayed_fmax_eV_A=observation['replayed_fmax_eV_A'],
            sampled_initial_peak_meV_fu=observation['sampled_forward_barrier_meV_fu'],
            observation_sha256=sha256(folder/'step_0000/observation.json'),
            new_DFT_calls_for_export=0, new_SCFs_in_this_initial_evaluation=7,
            terminal_path_claimed=False, old_holdout_labels_replaced=False, stationary_TS_certified=False)
        save(folder/'initial_receipt.json', report); results.append(report)
    files = {p.relative_to(output).as_posix(): sha256(p) for p in sorted(output.rglob('*')) if p.is_file()}
    save(output/'files.json', files)
    with tarfile.open(archive, 'x:gz') as tar:
        tar.add(output, arcname=output.name)
    receipt = dict(archive=str(archive), archive_sha256=sha256(archive), archive_bytes=archive.stat().st_size,
        exported_files=len(files), initial_evaluations=results, new_DFT_calls=0, scheduler_mutations=0)
    save(ROOT/'initial_collection_receipt.json', receipt)
    return receipt


if __name__ == '__main__':
    print(json.dumps(collect(), allow_nan=False))
