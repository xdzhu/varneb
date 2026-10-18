"""Inert exact-cache seed for one free-interior connection test; no DFT."""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import sys

REPOSITORY = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPOSITORY))
from ase.io import write
from scripts.audit_hfo2_static_replica import sha256

CASE = Path(__file__).resolve().parent
PARENT = CASE.with_name('preserving_pattern_branches_E272_20261018')
REMOTE_PARENT = '/public/home/iai806/abacus/agent-runs/20261018-varneb-preserving-branches-E272-r1'
ASSETS = ('/public/home/iai806/abacus/agent-runs/20260912-vcneb-hf/validation/hfo2_t_to_po/'
          'vcneb_n7_fire_distributed_resume27678004_v4/00')


def save(path, value):
    with path.open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(value, handle, indent=2, allow_nan=False); handle.write('\n')


def prepare(output):
    spec = importlib.util.spec_from_file_location('E273_original_connection', PARENT/'analyze_prefix.py')
    reader = importlib.util.module_from_spec(spec); spec.loader.exec_module(reader)
    report, images = reader.analyze(); _, original = reader.original_images()
    if (report['new_DFT_calls'] != 0 or report['replaced_image'] != 4
            or report['replaced_central_image']['ordinary_residual_pass']):
        raise ValueError('the native, not-yet-converged connection polygon is required')
    output = Path(output); output.mkdir(parents=True, exist_ok=False)
    records = [{k: row[k] for k in ('input_sha256', 'raw_log_sha256')}
               | {'directory': row['raw_source']} for row in original['raw_image_evaluations']]
    native = report['signed_branch']['arms']['J']
    record = Path(native['native_record_relative_path'])
    if record.parts[0] != 'evaluations' or record.name != 'result.json' or '..' in record.parts:
        raise ValueError('one exact originally evaluated J source is required')
    records[4] = dict(directory=f'{REMOTE_PARENT}/results/0/minus/{record.parent.as_posix()}/scf_000000',
        input_sha256=native['native_input_sha256'], raw_log_sha256=native['native_log_sha256'])
    save(output/'factory_parameters.json', dict(source_directory=ASSETS, seed_cache_records=records))
    write(output/'seed.traj', images)
    for name, image in (('initial.vasp', images[0]), ('final.vasp', images[-1]), ('substrate.vasp', images[0])):
        write(output/name, image, format='vasp', direct=True, vasp5=True)
    save(output/'connection_replay.json', report)
    files = {p.name: sha256(p) for p in sorted(output.iterdir()) if p.is_file()}
    manifest = dict(status='prepared_only_no_runtime_or_submission', source_archive_sha256=reader.ARCHIVE_SHA,
        source_observation_sha256=reader.OBSERVATION_SHA, branch='zero_minus_J', n_total_images=9,
        n_fixed_endpoints=2, free_interiors=list(range(1, 8)), retained_pattern_constraints_in_NEB=[],
        source_job_must_finish='28948425', original_labels_and_scores_unchanged=True,
        required_factory='examples.hfo2_fixed_input_factory:make_clamped_resume_cached_factory',
        fmax_target_eV_A=.10, maximum_updates=24, maximum_new_interior_SCFs=168,
        walltime='06:00:00', maximum_study_allocations=2, pressure_GPa=0., climb=False,
        new_DFT_calls=0, scheduler_mutations=0, execution_admitted=False, files_sha256=files,
        preparation_source_sha256=sha256(Path(__file__)), protocol_sha256=sha256(CASE/'protocol.md'))
    save(output/'manifest.json', manifest)
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); print(json.dumps(prepare(args.output)))
