"""Stage a separate B4-side pressure refinement from the audited plus descent.

The 1000 eV physical settings and 45.7 GPa pressure remain fixed. Only the
VASP force stopping criterion and maximum ionic step count change. This does
not alter the original 600 eV VCNEB path or the completed plus descent.
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

from ase.io import read

from scripts.prepare_gan_ts_basin_continuation import geometry_metrics
from scripts.prepare_gan_ts_newton_probe import same_geometry, sha256


def prepare(source_root: Path, source_audit_path: Path,
            structural_match_path: Path, work_root: Path) -> dict:
    if work_root.exists():
        raise FileExistsError(work_root)
    manifest_path = source_root / 'manifest.json'
    manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
    audit = json.loads(source_audit_path.read_text(encoding='utf-8'))
    match = json.loads(structural_match_path.read_text(encoding='utf-8'))
    if (manifest.get('purpose')
            != 'GaN_plus_basin_100step_uninterrupted_from_audited_seed_not_endpoint_certificate'
            or audit.get('case_kind') != 'uninterrupted_plus'
            or audit.get('source_sha256', {}).get('manifest') != sha256(manifest_path)
            or audit.get('first_ten_max_H_difference_from_pilot_eV', 1) > 1e-5
            or not audit.get('vasp_reports_ionic_convergence')
            or not audit.get('contcar_is_last_evaluated_geometry')
            or not audit.get('force_pass_0p02eV_per_A')
            or audit.get('stress_pass_1kbar')
            or match.get('endpoint') != 'B4'
            or match.get('candidate_spacegroup_at_0p01A') != 'P6_3mc'
            or match.get('candidate_coordination_GaN_2p4A') != [4, 4]
            or match.get('source_sha256', {}).get('candidate')
            != audit.get('source_sha256', {}).get('CONTCAR')):
        raise ValueError('source is not the audited B4 candidate needing pressure refinement')
    source = source_root / 'case'
    if (sha256(source / 'OUTCAR') != audit['source_sha256']['OUTCAR']
            or sha256(source / 'CONTCAR') != audit['source_sha256']['CONTCAR']
            or not same_geometry(read(source / 'CONTCAR', format='vasp'),
                                 read(source / 'OUTCAR'), tolerance=2e-5)):
        raise ValueError('B4 candidate geometry/output changed')
    geometry = geometry_metrics(read(source / 'CONTCAR', format='vasp'))
    incar = (source / 'INCAR').read_bytes()
    if (b'\r' in incar or incar.count(b'NSW = 100\n') != 1
            or incar.count(b'EDIFFG = -0.02\n') != 1
            or incar.count(b'PSTRESS = 457.0\n') != 1):
        raise ValueError('source VASP relaxation contract differs')
    updated = incar.replace(b'NSW = 100\n', b'NSW = 30\n')
    updated = updated.replace(b'EDIFFG = -0.02\n', b'EDIFFG = -0.005\n')
    destination = work_root / 'case'
    destination.mkdir(parents=True)
    shutil.copy2(source / 'CONTCAR', destination / 'POSCAR')
    for filename in ('KPOINTS', 'POTCAR'):
        shutil.copy2(source / filename, destination / filename)
    (destination / 'INCAR').write_bytes(updated)
    hashes = {filename: sha256(destination / filename)
              for filename in ('POSCAR', 'INCAR', 'KPOINTS', 'POTCAR')}
    (destination / 'sha256.inputs.json').write_text(
        json.dumps(hashes, indent=2) + '\n', encoding='utf-8')
    result = {
        'purpose': 'GaN_1000eV_B4_candidate_pressure_refinement_not_TS_certificate',
        'status': 'inputs_finalized_no_DFT', 'pressure_GPa': 45.7,
        'PSTRESS_kbar': 457.0, 'encut_eV': 1000,
        'changed': ['NSW 100 -> 30', 'EDIFFG -0.02 -> -0.005 eV/A'],
        'unchanged': ['POSCAR from last evaluated B4 candidate CONTCAR',
                      'KPOINTS Gamma 8x8x6', 'Ga_d+N POTCAR', 'ISYM -1',
                      'SYMPREC 1e-4', 'PSTRESS 457 kbar', 'IBRION 2', 'ISIF 3'],
        'initial_geometry_preflight': geometry,
        'input_sha256': hashes,
        'source_sha256': {'source_manifest': sha256(manifest_path),
                          'source_audit': sha256(source_audit_path),
                          'structural_match': sha256(structural_match_path),
                          'source_OUTCAR': sha256(source / 'OUTCAR'),
                          'source_CONTCAR': sha256(source / 'CONTCAR'),
                          'preparer': sha256(Path(__file__))},
        'limitations': [
            'A stricter force criterion may not guarantee the 1 kbar cell-stress target.',
            'Final raw E/forces/stress and a same-setting static point require separate audit.',
        ],
    }
    (work_root / 'manifest.json').write_text(
        json.dumps(result, indent=2) + '\n', encoding='utf-8')
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('source-root', 'source-audit', 'structural-match', 'work-root'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.source_root, args.source_audit,
                     args.structural_match, args.work_root)
    print(json.dumps({'status': result['status'], 'changed': result['changed']}))


if __name__ == '__main__':
    main()
