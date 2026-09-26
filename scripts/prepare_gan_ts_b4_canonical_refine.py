"""Canonicalize only an audited near-B4 basin endpoint for VASP restart.

VASP 6.3.2 rejects the otherwise-converged B4-side CONTCAR at its Bravais
classifier, before SCF. The official POSCAR guidance permits high-precision
symmetry reconstruction of a near-symmetric endpoint. This isolated diagnostic
does not symmetrize any VCNEB image or silently change SYMPREC/DFT settings.
The input perturbation and first evaluated enthalpy must be audited separately.
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read, write

from scripts.prepare_gan_ts_basin_continuation import geometry_metrics
from scripts.prepare_gan_ts_newton_probe import sha256
from scripts.compare_gan_basin_endpoint import spacegroup_symbol


def canonicalize_near_b4(candidate: Atoms, reference: Atoms) -> tuple[Atoms, dict]:
    """Project only a tightly screened B4 basin candidate onto its P6₃mc chart."""

    if (candidate.get_chemical_symbols() != ['Ga', 'Ga', 'N', 'N']
            or reference.get_chemical_symbols() != candidate.get_chemical_symbols()
            or spacegroup_symbol(candidate, 0.01) != 'P6_3mc'
            or spacegroup_symbol(reference, 0.01) != 'P6_3mc'):
        raise ValueError('source/reference are not mapped four-atom B4 endpoints')
    reference_cell = reference.cell.array
    target_cell = reference_cell.copy()
    lengths = candidate.cell.lengths()
    a_length = float((lengths[0] + lengths[1]) / 2)
    target_cell[0:2] *= a_length / np.linalg.norm(reference_cell[0])
    target_cell[2] *= float(lengths[2] / np.linalg.norm(reference_cell[2]))
    if (abs(np.linalg.norm(reference_cell[0]) - np.linalg.norm(reference_cell[1])) > 1e-7
            or abs(reference.cell.angles()[2] - 120) > 1e-6
            or np.max(np.abs(reference.cell.angles()[:2] - 90)) > 1e-6):
        raise ValueError('reference B4 cell is not an exact hexagonal chart')
    reference_scaled = reference.get_scaled_positions(wrap=False)
    candidate_scaled = candidate.get_scaled_positions(wrap=False)
    delta = candidate_scaled - reference_scaled
    delta -= np.rint(delta)
    target_scaled = reference_scaled.copy()
    target_scaled[:, :2] += np.mean(delta[:, :2], axis=0)
    for pair in ((0, 1), (2, 3)):
        target_scaled[list(pair), 2] += np.mean(delta[list(pair), 2])
    target_scaled %= 1
    canonical = candidate.copy()
    canonical.calc = None
    canonical.set_cell(target_cell, scale_atoms=False)
    canonical.set_scaled_positions(target_scaled)
    position_delta = candidate.positions - canonical.positions
    position_delta -= np.rint(position_delta @ np.linalg.inv(target_cell)) @ target_cell
    metrics = {
        'max_atomic_position_change_A': float(np.max(np.linalg.norm(position_delta, axis=1))),
        'max_cell_vector_change_A': float(np.max(np.linalg.norm(
            candidate.cell.array - target_cell, axis=1))),
        'relative_volume_change': float(canonical.get_volume() / candidate.get_volume() - 1),
        'candidate_spacegroup_at_0p01A': spacegroup_symbol(candidate, 0.01),
        'canonical_spacegroup_at_1e_minus_4A': spacegroup_symbol(canonical, 1e-4),
    }
    if (metrics['max_atomic_position_change_A'] > 0.005
            or metrics['max_cell_vector_change_A'] > 0.005
            or abs(metrics['relative_volume_change']) > 0.001
            or metrics['canonical_spacegroup_at_1e_minus_4A'] != 'P6_3mc'):
        raise ValueError(f'B4 endpoint canonicalization exceeds geometry gate: {metrics}')
    return canonical, metrics


def prepare(source_root: Path, source_audit_path: Path,
            structural_match_path: Path, original_path: Path,
            failed_refine_root: Path, work_root: Path) -> dict:
    if work_root.exists():
        raise FileExistsError(work_root)
    source_audit = json.loads(source_audit_path.read_text(encoding='utf-8'))
    match = json.loads(structural_match_path.read_text(encoding='utf-8'))
    failed_manifest = json.loads((failed_refine_root / 'manifest.json').read_text(encoding='utf-8'))
    failed_outcar = failed_refine_root / 'case/OUTCAR'
    if (source_audit.get('case_kind') != 'uninterrupted_plus'
            or not source_audit.get('vasp_reports_ionic_convergence')
            or not source_audit.get('contcar_is_last_evaluated_geometry')
            or match.get('endpoint') != 'B4'
            or failed_manifest.get('purpose')
            != 'GaN_1000eV_B4_candidate_pressure_refinement_not_TS_certificate'
            or failed_manifest.get('source_sha256', {}).get('source_audit')
            != sha256(source_audit_path)
            or not failed_outcar.is_file()
            or 'Inconsistent Bravais lattice types' not in failed_outcar.read_text(
                encoding='utf-8', errors='replace')
            or 'General timing and accounting informations' in failed_outcar.read_text(
                encoding='utf-8', errors='replace')):
        raise ValueError('not the audited B4 endpoint and parser-only failed refinement')
    source = source_root / 'case'
    if (sha256(source / 'OUTCAR') != source_audit['source_sha256']['OUTCAR']
            or sha256(source / 'CONTCAR') != source_audit['source_sha256']['CONTCAR']
            or sha256(source / 'CONTCAR') != match['source_sha256']['candidate']
            or sha256(failed_refine_root / 'case/POSCAR')
            != sha256(source / 'CONTCAR')):
        raise ValueError('B4 raw source or failed restart input changed')
    candidate = read(source / 'CONTCAR', format='vasp')
    reference = read(original_path, index=0)
    canonical, perturbation = canonicalize_near_b4(candidate, reference)
    destination = work_root / 'case'
    destination.mkdir(parents=True)
    write(destination / 'POSCAR', canonical, format='vasp', direct=True,
          sort=False, vasp5=True)
    for filename in ('INCAR', 'KPOINTS', 'POTCAR'):
        shutil.copy2(failed_refine_root / 'case' / filename, destination / filename)
    reread = read(destination / 'POSCAR', format='vasp')
    if (not np.allclose(reread.cell.array, canonical.cell.array, atol=1e-9, rtol=0)
            or not np.allclose(reread.positions, canonical.positions, atol=1e-8, rtol=0)
            or b'\r' in (destination / 'INCAR').read_bytes()
            or (destination / 'INCAR').read_bytes()
            != (failed_refine_root / 'case/INCAR').read_bytes()):
        raise ValueError('high-precision canonical B4 POSCAR or copied contract changed')
    hashes = {filename: sha256(destination / filename)
              for filename in ('POSCAR', 'INCAR', 'KPOINTS', 'POTCAR')}
    (destination / 'sha256.inputs.json').write_text(
        json.dumps(hashes, indent=2) + '\n', encoding='utf-8')
    result = {
        'purpose': 'GaN_1000eV_B4_canonical_endpoint_refine_after_VASP_Bravais_failure',
        'status': 'inputs_finalized_no_DFT',
        'pressure_GPa': 45.7, 'PSTRESS_kbar': 457.0,
        'DFT_settings_change': 'none',
        'geometry_change': 'explicit near-P6_3mc endpoint canonicalization, not a VCNEB image',
        'canonicalization_metrics': perturbation,
        'initial_geometry_preflight': geometry_metrics(reread),
        'input_sha256': hashes,
        'source_sha256': {'plus_audit': sha256(source_audit_path),
                          'B4_match': sha256(structural_match_path),
                          'plus_CONTCAR': sha256(source / 'CONTCAR'),
                          'failed_refine_manifest': sha256(failed_refine_root / 'manifest.json'),
                          'failed_refine_OUTCAR': sha256(failed_outcar),
                          'original_path': sha256(original_path),
                          'preparer': sha256(Path(__file__))},
        'limitations': [
            'The canonicalized POSCAR differs slightly from the last evaluated descent geometry.',
            'The first new enthalpy must be checked against the audited descent energy.',
            'A B4 basin confirmation is not by itself a full TS certificate.',
        ],
    }
    (work_root / 'manifest.json').write_text(json.dumps(result, indent=2) + '\n',
                                             encoding='utf-8')
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('source-root', 'source-audit', 'structural-match',
                 'original-path', 'failed-refine-root', 'work-root'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.source_root, args.source_audit, args.structural_match,
                     args.original_path, args.failed_refine_root, args.work_root)
    print(json.dumps({'status': result['status'],
                      'canonicalization_metrics': result['canonicalization_metrics']}))


if __name__ == '__main__':
    main()
