"""Inert native input preparation for the already locked E185 test points.

There is deliberately no submit/run mode: complete B3/released-curvature
evidence and final model forecasts must precede independent DFT labels.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
if (Path(__file__).resolve().parents[4] / 'vcneb').is_dir():
    sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np
from scripts import hfo2_G3_conditional_faces as faces
from examples.hfo2_fixed_input_factory import write_structure, geometry_roundtrip

study = faces.study
OLD = study.RUNS / '20261015-varneb-G3-faces-E161-r1/source'
RUNTIME_SHA = 'd07bb011c9e313afe23ca4e78207356ff98adfc43ca231bdd047e6d5d0f1be20'
FORECAST_SHA = 'aa77cdd7460ad809a5fbe06e80c684f7f73ce118835067c380ec219b96abd2d8'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def forecast(bundle):
    return faces.canary.pinned(bundle / 'forecast_v2.json', FORECAST_SHA)


def checked_start(data, point):
    """Use the literal frozen 40-vector; do not rebuild it by a projection."""
    x = np.asarray(point['starting_full_coordinates_A'])
    u = np.asarray(point['starting_internal_coordinates_A'])
    if (x.shape != (40,) or u.shape != (36,) or not np.isfinite(x).all()
            or not np.isfinite(u).all() or x[-1] != 0.
            or not np.allclose(data['basis'] @ u, x, atol=2e-12, rtol=0)
            or not np.allclose(data['axes'].T @ u, point['q_A'], atol=2e-12, rtol=0)
            or np.linalg.norm(u) > .5 or point['training_parent_id'] != 0
            or not faces.canary.geometry_guard(data['chart'], data['basis'], u)):
        raise ValueError('frozen start violates fixed Q/span/control/original geometry; no repair')
    atoms = data['chart'].to_atoms(x)
    if not np.array_equal(atoms.cell.array[:2], data['chart'].reference_atoms.cell.array[:2]):
        raise ValueError('frozen substrate changed')
    return atoms


def prepare(bundle):
    locked = forecast(bundle)
    expected = {(i*.05, j*.05) for i in (-1,0,1) for j in (-1,0,1) if (i,j) != (0,0)}
    channels = []
    if [r['channel'] for r in locked['channels']] != ['T','preserving']:
        raise ValueError('only the original two local centers are authorized')
    for row in locked['channels']:
        data = faces.load_inputs(row['channel'])
        if len(row['targets']) != 8 or {tuple(p['q_A']) for p in row['targets']} != expected:
            raise ValueError('exact eight reserved coordinates required')
        for point in row['targets']:
            checked_start(data, point)
            if point['independent_DFT_label'] is not None:
                raise ValueError('initial forecast must precede all independent labels')
        channels.append(dict(channel=row['channel'], points=[{
            k:p[k] for k in ('id','q_A','training_parent_id','starting_internal_coordinates_A',
                            'starting_full_coordinates_A')} for p in row['targets']]))
    return dict(schema_version=1, purpose='inert_initial_decks_not_DFT_admission',
        channels=channels, forecast_sha256=FORECAST_SHA,
        original_input_artifact_pins=faces.PINS, six_asset_sha256=study.CONTRACT,
        binary_sha256=study.BINARY_SHA, original_runtime_manifest_sha256=RUNTIME_SHA,
        helper_text_sha256=study.text_sha(Path(__file__)),
        future_attempt_limit_per_point=12, future_attempt_limit_total=192,
        independent_points_reserved_each=8, independent_points_computed_each=0,
        future_affine_evaluation='literal_frozen_full_start + B @ (u - frozen_internal_start)',
        future_single_start=True, future_gradient_target_eV_A=.03,
        maximum_internal_radius_A=.5, minimum_distance_A=1.6,
        no_test_to_test_continuation=True, production_test_admitted=False,
        physical_plus005_generated_or_read=False)


def native_preflight(config, bundle, workdir):
    if (study.ROOT.resolve() != OLD.resolve() or workdir.exists()
            or workdir.resolve() != bundle.resolve().parent / 'initial_decks'
            or not workdir.resolve().is_relative_to(study.RUNS)
            or any(os.environ.get(k) != '1' for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'))):
        raise ValueError('immutable HF runtime, single-thread replay and fresh dedicated inert deck root required')
    if config != json.loads(json.dumps(prepare(bundle))):
        raise ValueError('frozen input/geometry/forecast/helper configuration differs')
    files = faces.canary.pinned(OLD / 'runtime_files.json', RUNTIME_SHA)
    if len(files) != 380 or any(sha(OLD / n) != h for n,h in files.items()):
        raise ValueError('original 380-file runtime changed')
    for name, digest in study.CONTRACT.items():
        if sha(study.ASSETS / name) != digest:
            raise ValueError('original six physical assets changed')
    if sha(study.BINARY) != study.BINARY_SHA:
        raise ValueError('original ABACUS binary changed')
    data = {c:faces.native_preflight(faces.prepare(c), c) for c in ('T','preserving')}
    workdir.mkdir()
    decks = []
    for channel in config['channels']:
        for point in channel['points']:
            atoms = checked_start(data[channel['channel']], point)
            directory = workdir / channel['channel'] / f"point_{point['id']:03d}"
            directory.mkdir(parents=True)
            for name in study.CONTRACT:
                if name in ('INPUT','KPT'):
                    shutil.copyfile(study.ASSETS / name, directory / name)
                else:
                    (directory / name).symlink_to(study.ASSETS / name)
            write_structure(directory / 'STRU', atoms)
            geometry_roundtrip(directory / 'STRU', atoms)
            hashes = {name:sha(directory / name) for name in (*study.CONTRACT,'STRU')}
            if any(hashes[name] != digest for name,digest in study.CONTRACT.items()):
                raise ValueError('inert input serialization changed physical files')
            distances = atoms.get_all_distances(mic=True)
            np.fill_diagonal(distances,np.inf)
            row = dict(channel=channel['channel'], id=point['id'], q_A=point['q_A'],
                directory=str(directory), input_sha256=hashes,
                volume_A3=float(atoms.get_volume()), minimum_distance_A=float(distances.min()),
                native_structure_roundtrip_passed=True, DFT_started=False)
            study.exclusive_json(directory / 'inert_input_manifest.json', row)
            decks.append(row)
    queue = subprocess.run(['squeue','-h','-j','28867177,28868827','-o','%i|%T|%R'],
                           check=True,capture_output=True,text=True).stdout.strip()
    return dict(status='sixteen_inert_native_initial_decks_passed_not_production_admission',
        config_sha256=sha(bundle / 'config.json'), forecast_sha256=FORECAST_SHA,
        original_runtime_files_checked=380, old_native_origins_reparsed=2,
        decks=decks, existing_job_queue=queue, new_DFT_or_scheduler_writes=0,
        independent_points_computed=0, production_test_admitted=False,
        final_B3_released_curvature_and_candidate_forecast_gates_still_required=True,
        physical_plus005_generated_or_read=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('prepare','preflight'))
    parser.add_argument('--bundle', type=Path, required=True)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--config-sha256')
    parser.add_argument('--workdir', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.mode == 'prepare':
        study.exclusive_json(args.config, prepare(args.bundle))
        return
    if sha(args.config) != args.config_sha256:
        raise ValueError('exact frozen config bytes required')
    if args.workdir is None or args.output is None or args.output.exists():
        raise FileExistsError('fresh deck root and report required')
    result = native_preflight(json.loads(args.config.read_bytes()), args.bundle, args.workdir)
    study.exclusive_json(args.output, result)
    print(json.dumps(dict(inert_decks=len(result['decks']), new_DFT=0, production_admitted=False)))


if __name__ == '__main__':
    main()
