"""Audit raw VASP output from isolated canonical B4 endpoint refinement."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from scripts.audit_gan_ts_basin_followups import residual_stress_kbar
from scripts.audit_gan_ts_basin_pilot import energy_triples, ga_n_coordination
from scripts.prepare_gan_ts_newton_probe import same_geometry, sha256


def audit(work_root: Path, source_root: Path, source_audit_path: Path,
          structural_match_path: Path, failed_refine_root: Path,
          original_path: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    manifest_path = work_root / 'manifest.json'
    manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
    source_audit = json.loads(source_audit_path.read_text(encoding='utf-8'))
    match = json.loads(structural_match_path.read_text(encoding='utf-8'))
    source = source_root / 'case'
    failed_outcar = failed_refine_root / 'case/OUTCAR'
    expected_sources = {
        'plus_audit': source_audit_path,
        'B4_match': structural_match_path,
        'plus_CONTCAR': source / 'CONTCAR',
        'failed_refine_manifest': failed_refine_root / 'manifest.json',
        'failed_refine_OUTCAR': failed_outcar,
        'original_path': original_path,
        'preparer': Path(__file__).with_name('prepare_gan_ts_b4_canonical_refine.py'),
    }
    if (manifest.get('purpose')
            != 'GaN_1000eV_B4_canonical_endpoint_refine_after_VASP_Bravais_failure'
            or manifest.get('status') != 'inputs_finalized_no_DFT'
            or manifest.get('DFT_settings_change') != 'none'
            or manifest.get('PSTRESS_kbar') != 457.0
            or source_audit.get('case_kind') != 'uninterrupted_plus'
            or match.get('endpoint') != 'B4'
            or any(manifest.get('source_sha256', {}).get(key) != sha256(path)
                   for key, path in expected_sources.items())):
        raise ValueError('canonical B4 refinement provenance changed')
    case = work_root / 'case'
    if (any(sha256(case / name) != digest
            for name, digest in manifest['input_sha256'].items())
            or b'\r' in (case / 'INCAR').read_bytes()
            or (case / 'INCAR').read_bytes()
            != (failed_refine_root / 'case/INCAR').read_bytes()):
        raise ValueError('canonical B4 input contract changed')
    outcar = case / 'OUTCAR'
    raw = outcar.read_text(encoding='utf-8', errors='replace')
    if ('General timing and accounting informations' not in raw
            or 'PSTRESS=  457.0' not in raw
            or 'aborting loop because EDIFF is reached' not in raw):
        raise ValueError('canonical B4 refinement lacks complete converged VASP output')
    triples = energy_triples(raw)
    if not 1 <= len(triples) <= 30:
        raise ValueError('unexpected canonical B4 ionic energy record count')
    atoms = read(outcar)
    contcar = case / 'CONTCAR'
    candidate = read(contcar, format='vasp')
    if (len(atoms) != 4 or not np.isfinite(atoms.get_forces()).all()
            or not np.isfinite(atoms.get_stress(voigt=False)).all()
            or not np.isclose(triples[-1][2], 45.7 * GPa * atoms.get_volume(),
                              atol=2e-4, rtol=0)):
        raise ValueError('canonical B4 final E/forces/stress/PV invalid')
    first_delta = float(triples[0][1] - source_audit['last_enthalpy_eV_per_cell'])
    if abs(first_delta) > 0.002:
        raise ValueError(f'canonicalization changed initial enthalpy by {first_delta} eV')
    fmax = float(np.max(np.linalg.norm(atoms.get_forces(), axis=1)))
    stress = residual_stress_kbar(atoms.get_stress(voigt=False), 45.7)
    geometry_pass = same_geometry(atoms, candidate, tolerance=2e-5)
    ionic_pass = 'reached required accuracy - stopping structural energy minimisation' in raw
    force_pass = fmax <= 0.02
    stress_pass = stress <= 1.0
    result = {
        'status': ('GaN_B4_canonical_refine_raw_audited_force_stress_pass_not_TS_certificate'
                   if all((geometry_pass, ionic_pass, force_pass, stress_pass)) else
                   'GaN_B4_canonical_refine_raw_audited_review_required'),
        'endpoint_candidate': 'B4',
        'n_ionic_energy_records': len(triples),
        'initial_enthalpy_minus_raw_plus_endpoint_eV_per_cell': first_delta,
        'last_free_energy_eV_per_cell': triples[-1][0],
        'last_PV_eV_per_cell': triples[-1][2],
        'last_enthalpy_eV_per_cell': triples[-1][1],
        'last_volume_A3': float(atoms.get_volume()),
        'last_max_atomic_force_eV_per_A': fmax,
        'last_max_stress_residual_from_45p7GPa_kbar': stress,
        'last_coordination_GaN_2p4A': ga_n_coordination(atoms),
        'contcar_is_last_evaluated_geometry': geometry_pass,
        'vasp_reports_ionic_convergence': ionic_pass,
        'force_pass_0p02eV_per_A': force_pass,
        'stress_pass_1kbar': stress_pass,
        'source_sha256': {'manifest': sha256(manifest_path),
                          'source_audit': sha256(source_audit_path),
                          'structural_match': sha256(structural_match_path),
                          'OUTCAR': sha256(outcar), 'CONTCAR': sha256(contcar),
                          'auditor': sha256(Path(__file__))},
        'limitations': [
            'The canonicalization is an explicit sub-0.002 Å endpoint perturbation, not a path image change.',
            'B4 phase identity requires separate final structural comparison.',
            'A same-setting static E+PV check and other-side basin evidence remain separate gates.',
        ],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'status': result['status'], 'n_records': len(triples),
                      'H_eV': triples[-1][1], 'first_H_delta_eV': first_delta,
                      'fmax': fmax, 'stress_residual_kbar': stress}))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('work-root', 'source-root', 'source-audit', 'structural-match',
                 'failed-refine-root', 'original-path', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    audit(args.work_root, args.source_root, args.source_audit,
          args.structural_match, args.failed_refine_root,
          args.original_path, args.output)


if __name__ == '__main__':
    main()
