"""Run only the three preregistered visible endpoint features; no path SCFs."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys

import numpy as np
from ase.io import read

from examples.hfo2_fixed_input_factory import CONTRACT, make_factory, read_fixed_hfo2_stru, same_ordered_geometry
from scripts.audit_hfo2_clamped_canary import replay_directories
from scripts.prepare_hfo2_clamped_endpoints import geometry_guard
from scripts.prepare_hfo2_channel_work_probes import write_structure
from scripts.prepare_hfo2_gamma_and_seeds import geometry_roundtrip
from scripts.relax_clamped_ase_endpoint import load_seed, relax

ROOT = Path('/public/home/iai806/abacus/agent-runs/20261018-varneb-physical-E268-r1')
BIN = Path('/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus')
PHYSICS = Path('/public/home/iai806/abacus/agent-runs/20260912-vcneb-hf/validation/'
               'hfo2_t_to_po/vcneb_n7_fire_distributed_resume27678004_v4/00')
PHASES = (('PO_plus', 'PO_minus_T_preserving'), ('T',))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, obj):
    with Path(path).open('x', encoding='utf-8', newline='\n') as f:
        json.dump(obj, f, indent=2, allow_nan=False); f.write('\n')


def validate_bundle(root, manifest_sha):
    bundle = root/'bundle'
    path = bundle/'bundle_files.json'
    if sha(path) != manifest_sha:
        raise ValueError('externally pinned bundle changed')
    manifest = json.loads(path.read_bytes())
    for name, digest in manifest.items():
        relative = Path(name)
        if relative.is_absolute() or '..' in relative.parts or sha(bundle/relative) != digest:
            raise ValueError('bundle input drift')
    config = json.loads((bundle/'runtime.json').read_bytes())
    if sha(root/'source.tar') != config['source_archive_sha256']:
        raise ValueError('immutable source archive changed')
    for name, digest in config['source_files_sha256'].items():
        relative = Path(name)
        if relative.is_absolute() or '..' in relative.parts or sha(root/'source'/relative) != digest:
            raise ValueError('executed source changed')
    if any(sha(PHYSICS/name) != h for name, h in CONTRACT.items()):
        raise ValueError('six physics hashes changed')
    if sha(BIN) != '88ea7f91c9f3091ef5fb4f410ab88267b14ecfe7f35c86dfaee9bf9319cace7a':
        raise ValueError('ABACUS binary changed')
    if sha(bundle/'registration.json') != 'faebe7588579f5dca9253e672755d4db5c85f8a33165685dc5cdce3cef2899ff':
        raise ValueError('pre-endpoint registration changed')
    registration = json.loads((bundle/'registration.json').read_bytes())
    if (registration['physical_input_sha256'] != CONTRACT
            or registration['protocol_sha256'] != sha(bundle/'protocol.md')
            or registration['endpoint_shards'] != [list(p) for p in PHASES]):
        raise ValueError('physics, protocol or phases differ')
    return config


def preflight(root, pin):
    validate_bundle(root, pin)
    target = root/'preflight'
    target.mkdir(exist_ok=False)
    decks = []
    for phase in PHASES[0]+PHASES[1]:
        seed = root/'bundle/endpoint_seeds'/phase/'endpoint_seed.json'
        atoms, boundary, manifest = load_seed(seed)
        if (manifest['strain'] != .005 or manifest['phase_label'] != phase
                or manifest['pressure_gpa'] != 0 or manifest['physical_contract_sha256'] != CONTRACT
                or manifest['old_endpoint_results_reused'] or manifest['target_interior_generated']):
            raise ValueError('only predeclared target endpoint seeds allowed')
        deck = target/phase; deck.mkdir()
        for name in ('INPUT', 'KPT'): shutil.copyfile(PHYSICS/name, deck/name)
        write_structure(deck/'STRU', atoms); geometry_roundtrip(deck/'STRU', atoms)
        reread = read_fixed_hfo2_stru(deck/'STRU')
        boundary.validate_images([reread]); geometry_guard([reread])
        if not same_ordered_geometry(atoms, reread):
            raise ValueError('native roundtrip changed geometry')
        decks.append(dict(phase=phase, seed_manifest_sha256=sha(seed),
                          input_sha256={name: sha(deck/name) for name in ('INPUT', 'KPT', 'STRU')}))
    result = dict(status='three_native_endpoint_decks_ready', bundle_sha256=pin,
                  new_DFT_calls=0, target_interiors_generated=False, decks=decks)
    save(root/'native_preflight.json', result)
    return result


def audit_endpoint(root, phase):
    seed = root/'bundle/endpoint_seeds'/phase/'endpoint_seed.json'
    _, boundary, manifest = load_seed(seed)
    folder = root/'endpoints'/phase
    summary = json.loads((folder/'endpoint_relax_summary.json').read_bytes())
    calls = sorted((folder/'calculator/image_0000').glob('scf_*'))
    if (not calls or len(calls) != summary['optimizer_steps']+1 or len(calls) > 81
            or summary['seed_manifest_sha256'] != sha(seed) or summary['steps_requested'] != 80
            or summary['fmax_target_eV_per_A'] != .03 or summary['open_stress_target_kbar'] != 2.
            or summary['maxstep'] != .02 or summary['strain'] != .005
            or summary['physical_contract_sha256'] != CONTRACT):
        raise ValueError('endpoint recipe/count/source differs')
    rows = replay_directories(calls, boundary, full_physical_bytes=True)
    last = rows[-1]
    for key, value in [('potential_energy_eV', last['energy_eV_cell']),
                       ('max_atomic_force_eV_per_A', last['max_atomic_force_eV_A']),
                       ('open_traction_norm_kbar', last['open_traction_norm_kbar'])]:
        if not np.isclose(summary[key], value, atol=2e-8, rtol=0):
            raise ValueError('native terminal result differs')
    converged = last['max_atomic_force_eV_A'] < .03 and last['open_traction_norm_kbar'] < 2.
    if summary['converged'] != converged:
        raise ValueError('false endpoint convergence')
    if not same_ordered_geometry(read_fixed_hfo2_stru(calls[-1]/'STRU'), read(folder/'CONTCAR', format='vasp')):
        raise ValueError('terminal CONTCAR differs')
    result = dict(status='native_endpoint_audited_phase_review_required', phase=phase,
        strain=.005, endpoint_converged=converged, new_DFT_calls_for_analysis=0,
        new_SCF_calls=len(rows), real_MPI_ranks=32,
        raw_full_physical_bytes_checked_here=True, rows=rows,
        SCF_core_hours=sum(row['elapsed_seconds'] for row in rows)*32/3600,
        summary_sha256=sha(folder/'endpoint_relax_summary.json'), seed_manifest_sha256=sha(seed),
        cache=dict(directory=str(calls[-1]), input_sha256=last['input_sha256'],
                   raw_log_sha256=last['raw_log_sha256']),
        target_interiors_evaluated=False, phase_stability_certified=False)
    save(folder/'native_audit.json', result)
    return result


def run(root, pin, shard):
    validate_bundle(root, pin)
    if (os.environ.get('SLURM_JOB_PARTITION') != 'hfacnormal01'
            or os.environ.get('SLURM_JOB_NUM_NODES') != '1'
            or os.environ.get('SLURM_CPUS_PER_TASK') != '32'
            or os.environ.get('SLURM_NTASKS') != '1'
            or os.environ.get('SLURM_JOB_NAME') != f'hfo2-E268-endpoints-{shard}'):
        raise ValueError('one32CPU allocation on prescribed partition required')
    proof = json.loads((root/'native_preflight.json').read_bytes())
    if proof['bundle_sha256'] != pin or len(proof['decks']) != 3 or proof['new_DFT_calls'] != 0:
        raise ValueError('three-deck preflight required')
    accepted = json.loads((root/f'accepted_{shard}.json').read_bytes())
    if accepted['job_id'] != os.environ['SLURM_JOB_ID'] or accepted['bundle_sha256'] != pin:
        raise ValueError('only accepted original handle may execute')
    (root/'endpoints').mkdir(exist_ok=True)
    results = []
    for phase in PHASES[shard]:
        seed = root/'bundle/endpoint_seeds'/phase/'endpoint_seed.json'
        atoms, boundary, manifest = load_seed(seed)
        factory = make_factory(parameters={'source_directory': str(PHYSICS)},
                               command=f'mpirun -np 32 {BIN}')
        summary = relax(atoms, boundary, manifest, root/'endpoints'/phase,
                        factory=factory, fmax=.03, stress_kbar=2., steps=80,
                        maxstep=.02, seed_manifest_sha256=sha(seed))
        audit = audit_endpoint(root, phase)
        results.append(dict(phase=phase, converged=summary['converged'],
                            new_SCF_calls=audit['new_SCF_calls'], SCF_core_hours=audit['SCF_core_hours']))
        print(json.dumps(results[-1]), flush=True)
    save(root/f'endpoint_result_{shard}.json', dict(job_id=os.environ['SLURM_JOB_ID'],
        bundle_sha256=pin, results=results, target_interiors_evaluated=False))
    return results


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=('preflight', 'run'))
    p.add_argument('--bundle-sha256', required=True)
    p.add_argument('--shard', type=int, choices=(0, 1))
    a = p.parse_args()
    if a.mode == 'preflight': result = preflight(ROOT, a.bundle_sha256)
    else: result = run(ROOT, a.bundle_sha256, a.shard)
    print(json.dumps(result if a.mode == 'preflight' else {'completed_phases': len(result)}))
