"""Combine audited 1000 eV GaN local-index-one and bilateral basin evidence.

This reports pressure enthalpy barriers per GaN formula unit only after both
same-setting basin statics, phase checks and the two-step joint Hessian pass.
It does not retrofit these 1000 eV energies onto the original 600 eV VCNEB.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from scripts.prepare_gan_ts_newton_probe import sha256


FILES = {
    'hessian_0p02': 'refined_hessian_0p02_audit_v1/audit.json',
    'hessian_0p01': 'refined_hessian_0p01_audit_v1/audit.json',
    'hessian_comparison': 'refined_hessian_0p02_0p01_comparison_v2.json',
    'path_mode_overlap_0p02': 'refined_hessian_0p02_mode_path_overlap_v2.json',
    'path_mode_overlap_0p01': 'refined_hessian_0p01_mode_path_overlap_v1.json',
    'signed_statics': 'signed_descent_0p10_audit_v1.json',
    'plus_descent': 'basin_plus_uninterrupted_100step_v1_audit.json',
    'minus_descent': 'basin_continuation_40step_v1_minus_audit_v1.json',
    'B4_refine': 'basin_B4_canonical_refine_v1_audit.json',
    'B1_refine': 'basin_B1_pressure_refine_v1_audit.json',
    'B4_match': 'basin_B4_canonical_refine_v1_structural_match.json',
    'B1_match': 'basin_B1_pressure_refine_v1_structural_match.json',
    'B4_static': 'basin_B4_static_1000eV_v1_audit.json',
    'B1_static': 'basin_B1_static_1000eV_v1_audit.json',
    'B4_static_manifest': 'basin_B4_static_1000eV_v1/manifest.json',
    'B1_static_manifest': 'basin_B1_static_1000eV_v1/manifest.json',
}


def enthalpy_barriers_per_formula_unit(
    center_eV_per_cell: float,
    initial_eV_per_cell: float,
    final_eV_per_cell: float,
    formula_units_per_cell: int,
) -> tuple[float, float, float]:
    """Return forward/reverse activation enthalpies and final−initial enthalpy.

    All three inputs must already be *enthalpies* at the same pressure and
    electronic-structure settings. In particular, this function never adds PV.
    """
    if (formula_units_per_cell <= 0
            or any(not math.isfinite(value) for value in (
                center_eV_per_cell, initial_eV_per_cell, final_eV_per_cell))
            or center_eV_per_cell <= max(initial_eV_per_cell, final_eV_per_cell)):
        raise ValueError('invalid same-setting enthalpies or formula-unit count')
    scale = formula_units_per_cell
    return (
        (center_eV_per_cell - initial_eV_per_cell) / scale,
        (center_eV_per_cell - final_eV_per_cell) / scale,
        (final_eV_per_cell - initial_eV_per_cell) / scale,
    )


def summarize(root: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    paths = {key: root / relative for key, relative in FILES.items()}
    records = {key: json.loads(path.read_text(encoding='utf-8'))
               for key, path in paths.items()}
    a02, a01 = records['hessian_0p02'], records['hessian_0p01']
    comparison = records['hessian_comparison']
    signed = records['signed_statics']
    if (comparison.get('kind') != 'gan'
            or comparison.get('step_sizes') != [0.02, 0.01]
            or comparison.get('source_sha256', {}).get('first_audit') != sha256(paths['hessian_0p02'])
            or comparison.get('source_sha256', {}).get('second_audit') != sha256(paths['hessian_0p01'])
            or any(counts != [1, 1] for counts in comparison['negative_counts_by_cutoff'].values())
            or comparison['archived_mode_comparison']['lowest_mode_absolute_overlap'] < 0.999
            or a01['source_sha256']['center_OUTCAR'] != a02['source_sha256']['center_OUTCAR']
            or abs(a01['center_enthalpy_eV_per_cell']
                   - a02['center_enthalpy_eV_per_cell']) > 1e-9
            or max(a01['center_gradient_translation_free_eV_per_A'],
                   a02['center_gradient_translation_free_eV_per_A']) > 0.02
            or signed.get('source_sha256', {}).get('hessian_audit') != sha256(paths['hessian_0p01'])
            or signed.get('source_sha256', {}).get('step_comparison')
            != sha256(paths['hessian_comparison'])
            or not signed.get('both_sides_downhill')
            or signed.get('n_completed_raw_VASP_statics') != 2):
        raise ValueError('joint-index-one center or signed descent evidence incomplete')
    for key in ('path_mode_overlap_0p02', 'path_mode_overlap_0p01'):
        record = records[key]
        if (record.get('pressure_GPa') != 45.7
                or record.get('peak_image_index') != 15
                or record.get('absolute_unstable_mode_path_tangent_overlap', 0) < 0.99):
            raise ValueError(f'path tangent/negative mode evidence invalid: {key}')
    plus, minus = records['plus_descent'], records['minus_descent']
    if (plus.get('case_kind') != 'uninterrupted_plus'
            or minus.get('case_kind') != 'continued_minus'
            or plus.get('first_ten_max_H_difference_from_pilot_eV', 1) > 1e-5
            or not all(record.get('vasp_reports_ionic_convergence')
                       and record.get('contcar_is_last_evaluated_geometry')
                       and record.get('force_pass_0p02eV_per_A')
                       for record in (plus, minus))
            or plus['last_enthalpy_eV_per_cell'] >= signed['center_enthalpy_eV_per_cell']
            or minus['last_enthalpy_eV_per_cell'] >= signed['center_enthalpy_eV_per_cell']):
        raise ValueError('both raw negative-mode descents are not verified')
    for phase, descent in (('B4', plus), ('B1', minus)):
        refine = records[f'{phase}_refine']
        match = records[f'{phase}_match']
        static = records[f'{phase}_static']
        static_manifest = records[f'{phase}_static_manifest']
        if (match.get('endpoint') != phase
                or match.get('candidate_spacegroup_at_0p01A')
                != ('P6_3mc' if phase == 'B4' else 'Fm-3m')
                or match.get('candidate_coordination_GaN_2p4A')
                != ([4, 4] if phase == 'B4' else [6, 6])
                or match.get('relative_displacement_rms_A', 1) > 0.01
                or refine.get('last_max_atomic_force_eV_per_A', 1) > 0.02
                or refine.get('last_max_stress_residual_from_45p7GPa_kbar', 2) > 1
                or not refine.get('force_pass_0p02eV_per_A')
                or not refine.get('stress_pass_1kbar')
                or not refine.get('vasp_reports_ionic_convergence')
                or not refine.get('contcar_is_last_evaluated_geometry')
                or static.get('status')
                != 'GaN_basin_static_1000eV_raw_audited_force_stress_and_H_match'
                or static.get('endpoint_candidate') != phase
                or not static.get('H_agrees_with_relax_within_2meV')
                or static_manifest.get('endpoint') != phase
                or static_manifest.get('pressure_for_postprocessing_GPa') != 45.7
                or static_manifest.get('source_sha256', {}).get('relax_audit')
                != sha256(paths[f'{phase}_refine'])
                or static_manifest.get('source_sha256', {}).get('structural_match')
                != sha256(paths[f'{phase}_match'])):
            raise ValueError(f'{phase} basin static, phase or force/stress gate failed')
        if phase == 'B4' and (
                abs(refine.get('initial_enthalpy_minus_raw_plus_endpoint_eV_per_cell', 1)) > 0.002
                or refine['last_enthalpy_eV_per_cell'] >= descent['last_enthalpy_eV_per_cell']):
            raise ValueError('B4 explicit canonicalization lacks small-enthalpy link')
    input_hashes = [records[f'{phase}_static_manifest']['input_sha256']
                    for phase in ('B4', 'B1')]
    for key in ('INCAR', 'KPOINTS', 'POTCAR'):
        if (input_hashes[0][key] != input_hashes[1][key]
                or input_hashes[0][key] != signed['cases'][0]['input_sha256'][key]):
            raise ValueError(f'1000 eV center/basin DFT contract differs: {key}')
    H_ts = float(a01['center_enthalpy_eV_per_cell'])
    H_b4 = float(records['B4_static']['H_E0_plus_PV_eV_per_cell'])
    H_b1 = float(records['B1_static']['H_E0_plus_PV_eV_per_cell'])
    forward, reverse, relative = enthalpy_barriers_per_formula_unit(
        H_ts, H_b4, H_b1, formula_units_per_cell=2)
    result = {
        'status': 'GaN_1000eV_near_stationary_joint_index_one_two_basin_links_supported',
        'pressure_GPa': 45.7,
        'calculator_contract': 'VASP 6.3.2 PBE/Ga_d+N PAW ENCUT1000 Gamma8x8x6 ISYM-1 SYMPREC1e-4',
        'cell_composition': 'Ga2N2 = 2 GaN formula units',
        'center_translation_free_gradient_eV_per_A': a01['center_gradient_translation_free_eV_per_A'],
        'joint_hessian_lowest_two_eV_per_A2': comparison['lowest_eigenvalues'],
        'joint_hessian_next_two_eV_per_A2': comparison['next_eigenvalues'],
        'joint_negative_mode_absolute_step_overlap': comparison['archived_mode_comparison']['lowest_mode_absolute_overlap'],
        'negative_mode_path_tangent_absolute_overlap': records['path_mode_overlap_0p01']['absolute_unstable_mode_path_tangent_overlap'],
        'H_TS_candidate_eV_per_Ga2N2_cell': H_ts,
        'H_B4_eV_per_Ga2N2_cell': H_b4,
        'H_B1_eV_per_Ga2N2_cell': H_b1,
        'forward_B4_to_B1_barrier_eV_per_GaN': forward,
        'reverse_B1_to_B4_barrier_eV_per_GaN': reverse,
        'B1_minus_B4_enthalpy_eV_per_GaN': relative,
        'B4_canonicalization_initial_H_jump_eV_per_cell': records['B4_refine']['initial_enthalpy_minus_raw_plus_endpoint_eV_per_cell'],
        'source_sha256': {key: sha256(path) for key, path in paths.items()},
        'auditor_sha256': sha256(Path(__file__)),
        'limitations': [
            'These are same-setting 1000 eV diagnostic enthalpy barriers, not substitutions for the original 600 eV path energies.',
            'The B4 endpoint restart used an explicitly bounded 0.00141 Å canonicalization after a parser-only failure.',
            'Local index one and two basin links do not prove that this is the global minimum-enthalpy path or finite-temperature free-energy saddle.',
            'The 18-dimensional atomic/strain Hessian uses a declared cell-coordinate scale; fixed-cell Gamma phonons alone are insufficient.',
        ],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({key: result[key] for key in (
        'status', 'forward_B4_to_B1_barrier_eV_per_GaN',
        'reverse_B1_to_B4_barrier_eV_per_GaN',
    )}))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    summarize(args.root, args.output)


if __name__ == '__main__':
    main()
