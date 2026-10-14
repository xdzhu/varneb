"""Finite same-input HfO2 training probes at two non-certified bottlenecks.

Preparation is zero DFT. Run mode is restricted to one fresh HF/32-MPI
allocation per declared centre, with exactly148 reserved static attempts.
No centre is a TS merely because its parent band passed ordinary NEB.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.stress import voigt_6_to_full_3x3_stress

from examples.hfo2_fixed_input_factory import CONTRACT, make_factory, read_fixed_hfo2_stru, same_ordered_geometry
from scripts.analyze_hfo2_clamped_reference import terminal_images
from scripts.audit_hfo2_static_replica import audited_results, sha256
from vcneb import (BiaxialClampedModeCoordinates, CalculatorModeEvaluator,
    assemble_joint_directional_curvature, clamped_plane_vcneb_boundary, joint_curvature_probes)

ROOT = Path(__file__).resolve().parents[1]
RUNS = Path('/public/home/iai806/abacus/agent-runs')
ASSETS = RUNS/'20260912-vcneb-hf/validation/hfo2_t_to_po/vcneb_n7_fire_distributed_resume27678004_v4/00'
BINARY = Path('/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus')
BINARY_SHA = '88ea7f91c9f3091ef5fb4f410ab88267b14ecfe7f35c86dfaee9bf9319cace7a'
COMMAND = 'mpirun -np 32 '+str(BINARY)
SCALE = 5.12968067458423
STEPS = (.01, .02)
SOURCES = {
    'T': ('clamped_T_terminal_E079_20261011',
        'c5653b4759b202199cb75739f99038aaa52b3ccc220ef226efbacc980e5d716e', '28692775'),
    'preserving': ('clamped_held_flips_E100_20261011/zero_strain_preserving',
        '46145da14f24b2b8bd46aa4656f0982650d6b2dc661282f58d809aab8a6d3202', '28727382'),
}
CODE = ('scripts/hfo2_G3_local_curvature.py', 'vcneb/biaxial_mode_coordinates.py',
    'vcneb/biaxial_curvature.py', 'vcneb/mode_evaluator.py', 'vcneb/joint_stencil.py',
    'vcneb/epitaxial_boundary.py', 'vcneb/active_curvature.py', 'vcneb/joint_curvature.py',
    'vcneb/core.py', 'examples/hfo2_fixed_input_factory.py', 'scripts/audit_hfo2_static_replica.py')
CODE += ('vcneb/clamped_mode_coordinates.py', 'vcneb/abacus.py',
    'scripts/prepare_hfo2_channel_work_probes.py', 'scripts/prepare_hfo2_gamma_and_seeds.py')


def text_sha(path):
    import hashlib
    return hashlib.sha256(path.read_text(encoding='utf-8').encode('utf-8')).hexdigest()


def exclusive_json(path, value):
    with Path(path).open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(value, stream, indent=2, allow_nan=False); stream.write('\n')


def chart_for(record):
    ref = record['reference']
    atoms = Atoms(ref['symbols'], cell=ref['cell_A'], positions=ref['positions_A'], pbc=ref['pbc'])
    boundary = clamped_plane_vcneb_boundary(12, np.asarray(record['substrate_reference_cell_A']), allow_tilt=True)
    return BiaxialClampedModeCoordinates(atoms, boundary=boundary,
        cell_scale_A=SCALE, anchor_strain=0., strain_scale_A=SCALE)


def probes_for(record):
    chart = chart_for(record)
    basis = np.asarray(record['directions'], dtype=float)
    if (chart.size != 40 or basis.shape != (40, 37)
            or not np.allclose(basis.T@basis, np.eye(37), atol=1e-12, rtol=0)
            or not np.array_equal(basis[-1], np.r_[np.zeros(36), 1.])
            or not np.allclose(chart.rigid_translation_directions().T@basis, 0., atol=1e-12, rtol=0)):
        raise ValueError('full translation-free36internal plus one separate control direction required')
    probes = []
    for h in STEPS:
        # Refuse a reserved control before generating its geometry.
        if any(abs(sign*h/SCALE-.005) < 1e-12 for sign in (1., -1.)):
            raise ValueError('reserved physical strain is sealed')
        for point in joint_curvature_probes(chart, basis, step_A=h):
            distances = point.atoms.get_all_distances(mic=True)
            np.fill_diagonal(distances, np.inf)
            if distances.min() < 1.6 or np.linalg.cond(point.atoms.cell.array) > 1e12:
                raise ValueError('unsafe probe geometry; do not repair or call DFT')
            probes.append(point)
    if len(probes) != 148:
        raise ValueError('exact finite stencil required')
    return chart, probes


def prepare(output):
    from scripts.audit_hfo2_clamped_network import audit
    from scripts.audit_hfo2_manuscript_work_cohort import BASE
    if output.exists():
        raise FileExistsError('fresh immutable preparation target required')
    network = audit()  # requires all eight native ordinary terminals
    raw = network['response']['sampled_minima']['0.0']
    if raw['sampled_switching_candidate'] != 'preserving' or raw['sampled_escape_candidate'] != 'T':
        raise ValueError('observed training selection differs; review, never silently replace centres')
    records = {}
    for name, (case, pin, job) in SOURCES.items():
        images, obs, residual, receipt = terminal_images(BASE/case, pin)
        peak = int(np.argmax([a.get_potential_energy() for a in images]))
        if peak != 4 or receipt['job_id'] != job:
            raise ValueError('registered source/peak differs')
        center = images[peak]
        record = dict(case=case, terminal_audit_sha256=pin, source_job_id=job, peak_image=peak,
            ordinary_fmax_eV_A=residual['fmax_eV_A'], raw_source=obs['raw_image_evaluations'][peak],
            reference=dict(symbols=center.get_chemical_symbols(), cell_A=center.cell.array.tolist(),
                positions_A=center.positions.tolist(), pbc=center.pbc.tolist()),
            substrate_reference_cell_A=obs['mechanical_boundary']['reference_cell_A'])
        chart = chart_for(record)
        record['directions'] = np.column_stack((chart.internal_relaxation_basis(), chart.controlled_direction())).tolist()
        # Ordered neighbours already use the continuous production lift; no MIC.
        tangent = chart.from_atoms(images[5])-chart.from_atoms(images[3])
        if not np.isfinite(tangent).all() or np.linalg.norm(tangent) <= 0 or abs(tangent[-1]) > 1e-12:
            raise ValueError('invalid fixed-strain native tangent')
        record['ordered_local_tangent'] = (tangent/np.linalg.norm(tangent)).tolist()
        _, probes = probes_for(record)
        record['probe_coordinates_A'] = [p.delta_A.tolist() for p in probes]
        records[name] = record
    manifest = dict(schema_version=1, purpose='two_training_candidate_full_local_controlled_curvatures_not_TS',
        ordinary_G2_paths=8, physical_inputs_sha256=dict(CONTRACT), binary_sha256=BINARY_SHA,
        normalized_source_text_sha256={p:text_sha(ROOT/p) for p in CODE},
        pressure_GPa=0., anchor_strain=0., cell_scale_A=SCALE, strain_scale_A=SCALE,
        step_amplitudes_A=list(STEPS), planned_raw_SCF_attempts_per_center=148,
        total_planned_raw_SCF_attempts=296, maximum_concurrent_allocations=2,
        MPI_ranks_per_allocation=32, partition='hfacnormal01', allocation_wall_hours=12,
        selection='lowest sampled zero-strain switching/escape candidates; exploratory, not certified physical ordering',
        barrier_error_bounds_measured=False, final_predictive_edge_selection_locked=False,
        unknown_error_policy='No finite error interval is fabricated. Final predictive selection must still apply the original uncertainty/tie rule; no third curvature centre is authorized if representatives change.',
        new_production_chains=0, new_material_or_orientation=False, climb=False,
        physical_holdout_generated_or_read=False, conditional_faces_selected=False,
        grid_nodes_used=0, independent_face_test_points_used=0,
        curvature_budget_scope='two allowed G3 centres;296local derivative SCFs, separately charged from the later two81-node/eight-test conditional faces',
        reference_policy='full native gradient; no source ASR, author-LDA curvature, fitted mother phase or deleted Cmma direction',
        centers=records)
    output.parent.mkdir(parents=True, exist_ok=True)
    exclusive_json(output, manifest)
    return manifest


def validate(manifest, channel):
    expected = dict(schema_version=1, pressure_GPa=0., anchor_strain=0., cell_scale_A=SCALE,
        strain_scale_A=SCALE, step_amplitudes_A=list(STEPS), planned_raw_SCF_attempts_per_center=148,
        total_planned_raw_SCF_attempts=296, MPI_ranks_per_allocation=32, partition='hfacnormal01',
        maximum_concurrent_allocations=2, allocation_wall_hours=12, new_production_chains=0,
        physical_holdout_generated_or_read=False, climb=False, ordinary_G2_paths=8,
        barrier_error_bounds_measured=False, final_predictive_edge_selection_locked=False,
        conditional_faces_selected=False, grid_nodes_used=0, independent_face_test_points_used=0)
    if (set(manifest['centers']) != set(SOURCES) or channel not in SOURCES
            or any(type(manifest.get(k)) is not type(v) or manifest.get(k) != v for k,v in expected.items())
            or manifest.get('physical_inputs_sha256') != CONTRACT or manifest.get('binary_sha256') != BINARY_SHA
            or manifest.get('normalized_source_text_sha256') != {p:text_sha(ROOT/p) for p in CODE}):
        raise ValueError('registered numerical/physical/code/stencil contract differs')
    record = manifest['centers'][channel]
    case, pin, job = SOURCES[channel]
    if (record['case'], record['terminal_audit_sha256'], record['source_job_id'], record['peak_image']) != (case,pin,job,4):
        raise ValueError('centre/source identity differs')
    chart, probes = probes_for(record)
    if not np.array_equal([p.delta_A for p in probes], record['probe_coordinates_A']):
        raise ValueError('prepared stencil coordinates differ')
    return chart, probes


def native_preflight(manifest, channel):
    chart, probes = validate(manifest, channel)
    if sha256(BINARY) != BINARY_SHA or any(sha256(ASSETS/n) != h for n,h in CONTRACT.items()):
        raise ValueError('existing executable or original six physical assets changed')
    row = manifest['centers'][channel]['raw_source']; source = Path(row['raw_source'])
    if (sha256(source/'call_audit.json') != row['audit_sha256']
            or sha256(source/'OUT.ABACUS/running_scf.log') != row['raw_log_sha256']
            or any(sha256(source/n) != h for n,h in row['input_sha256'].items())
            or not same_ordered_geometry(chart.reference_atoms, read_fixed_hfo2_stru(source/'STRU'))):
        raise ValueError('native centre provenance/geometry changed')
    raw = audited_results(source)
    if (raw['energy'] != row['energy_eV_cell']
            or not np.array_equal(raw['forces'], row['forces_eV_A'])
            or not np.array_equal(raw['stress'], row['stress_ASE_voigt_eV_A3'])):
        raise ValueError('native EFS changed')
    center = chart.reference_atoms
    center.calc = SinglePointCalculator(center, **raw)
    energy, gradient = chart.energy_and_gradient(np.zeros(40), center)
    return chart, probes, energy, gradient


def result_validator(index, atoms, directory):
    rawdirs = list(directory.glob('scf_*'))
    if len(rawdirs) != 1 or rawdirs[0].name != 'scf_000000':
        raise ValueError('exactly one fresh native SCF per static probe required')
    source = rawdirs[0]; raw = audited_results(source)
    call = json.loads((source/'call_audit.json').read_text())
    elapsed = call.get('elapsed_seconds')
    if isinstance(elapsed, bool) or not isinstance(elapsed, (float, int)) or not np.isfinite(elapsed) or elapsed <= 0:
        raise ValueError('finite positive native call time required')
    if (any(sha256(source/n) != h or call['input_sha256'].get(n) != h for n,h in CONTRACT.items())
            or sha256(source/'STRU') != call['input_sha256']['STRU']
            or sha256(source/'OUT.ABACUS/running_scf.log') != call['raw_log_sha256']
            or not same_ordered_geometry(atoms, read_fixed_hfo2_stru(source/'STRU'))
            or raw['energy'] != atoms.get_potential_energy()
            or not np.array_equal(raw['forces'], atoms.get_forces())
            or not np.array_equal(voigt_6_to_full_3x3_stress(raw['stress']), atoms.get_stress(voigt=False))):
        raise ValueError('fresh native raw EFS/input/geometry audit differs')
    return dict(raw_directory=str(source), raw_log_sha256=call['raw_log_sha256'],
        call_audit_sha256=sha256(source/'call_audit.json'), input_sha256=call['input_sha256'],
        genuine_MPI_ranks=32, raw_SCFs=1, elapsed_seconds=call['elapsed_seconds'])


def run(manifest, channel, workdir):
    if (os.environ.get('RUN_DFT') != '1' or os.environ.get('SLURM_JOB_PARTITION') != 'hfacnormal01'
            or os.environ.get('SLURM_JOB_NUM_NODES') != '1' or os.environ.get('SLURM_NTASKS') != '1'
            or os.environ.get('SLURM_CPUS_PER_TASK') != '32'
            or any(os.environ.get(n) != '1' for n in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'))):
        raise ValueError('one authorized hfacnormal01/32CPU allocation and one thread per MPI required')
    if not workdir.resolve().is_relative_to(RUNS) or workdir.exists():
        raise ValueError('fresh HF agent-runs namespace required')
    chart, probes, energy, gradient = native_preflight(manifest, channel)
    workdir.mkdir(parents=True)
    exclusive_json(workdir/'center.json', dict(cached_energy_eV_cell=energy,
        cached_full_gradient_eV_A=gradient.tolist(), new_center_SCFs=0,
        not_stationary_or_TS_certified=True, job_id=os.environ['SLURM_JOB_ID']))
    evaluate = CalculatorModeEvaluator(chart,
        make_factory(parameters={'source_directory':str(ASSETS)}, command=COMMAND),
        workdir=workdir/'evaluations', calculator_id='original-six-HfO2-PBE-100Ry-full10auDZP-Gamma2x2x2',
        minimum_distance_A=1.6, result_validator=result_validator,
        validation_id='E150-native-six-input-DSIZE32-EFS-one-SCF-v1', max_evaluation_attempts=148)
    energies, gradients = [], []
    for index, probe in enumerate(probes):
        e, g = evaluate(probe.delta_A); energies.append(e); gradients.append(g)
        print(json.dumps(dict(completed=index+1, total=148, step_A=probe.step_A,
            direction=probe.direction_index, sign=probe.sign, energy_eV_cell=e)), flush=True)
    basis = np.asarray(manifest['centers'][channel]['directions']); blocks = []
    for j,h in enumerate(STEPS):
        g = np.asarray(gradients[j*74:(j+1)*74])
        block = assemble_joint_directional_curvature(basis,g[::2],g[1::2],step_A=h)
        blocks.append(dict(step_A=h, raw_projected_eV_A2=block.raw_projected.tolist(),
            full_hessian_action_eV_A2=block.full_hessian_action.tolist(),
            reciprocity_relative_defect=block.reciprocity_relative_defect,
            symmetrization_operator_change_eV_A2=block.symmetrization_operator_change_eV_A2,
            fixed_strain_internal_eigenvalues_eV_A2=np.linalg.eigvalsh(block.symmetric_projected[:36,:36]).tolist()))
    completed = sorted((workdir/'evaluations').glob('eval-*/result.json'))
    validations = [json.loads(p.read_text())['validation'] for p in completed]
    if len(completed) != 148 or evaluate.n_evaluation_attempts != 148 or evaluate.n_new_evaluations != 148:
        raise ValueError('finite all-probe native completion differs')
    summary = dict(status='full_local_training_curvature_not_TS_or_material_prediction_certificate',
        source_channel=channel, native_probe_SCFs=148, new_production_chains=0,
        cell_scale_A=SCALE, strain_scale_A=SCALE, raw_center_gradient_eV_A=gradient.tolist(),
        energies_eV_cell=energies, full_gradients_eV_A=np.asarray(gradients).tolist(), blocks=blocks,
        sum_native_call_seconds=sum(v['elapsed_seconds'] for v in validations),
        call_time_proxy_core_hours=sum(v['elapsed_seconds'] for v in validations)*32/3600,
        allocation_core_hours=None, preparation_and_G2_costs_additional=True,
        fixed_strain_index_excludes_control=True, source_ASR_or_force_repair=False,
        physical_holdout_generated_or_read=False, conditional_surface_or_predictions_completed=False)
    exclusive_json(workdir/'curvature.json',summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=('prepare','preflight','run'))
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--manifest-sha256')
    parser.add_argument('--channel',choices=tuple(SOURCES))
    parser.add_argument('--workdir',type=Path)
    args=parser.parse_args()
    if args.mode == 'prepare':
        m=prepare(args.manifest); print(json.dumps(dict(centers=list(m['centers']),static_attempts=296)))
        return
    if sha256(args.manifest) != args.manifest_sha256 or args.channel is None:
        raise ValueError('exact frozen manifest digest and registered channel required')
    m=json.loads(args.manifest.read_text())
    if args.mode == 'preflight':
        native_preflight(m,args.channel); print('native preflight passed; zero new DFT')
    else:
        if args.workdir is None: raise ValueError('fresh workdir required')
        run(m,args.channel,args.workdir)


if __name__ == '__main__':
    main()
