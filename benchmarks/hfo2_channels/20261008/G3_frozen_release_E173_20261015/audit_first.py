"""One read-only native check of the first completed new matched-arm SCF.

Run with original E161 source on PYTHONPATH. Later unfinished trials are
pending, not failed. This never calls a calculator, resumes or submits.
"""
from datetime import datetime,timezone,timedelta
import hashlib,importlib.util,json,subprocess,sys
from pathlib import Path
sys.dont_write_bytecode=True
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator

RUNS=Path('/public/home/iai806/abacus/agent-runs')
ROOT=RUNS/'20261015-varneb-G3-feasible-execution-E171-r1'
OUT=RUNS/'20261015-varneb-G3-first-feasible-native-E173-r1'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()


def recorded_evaluation_directory(raw_directory,expected_evaluation,allowed_root):
    """Preserve recorded spelling only after strict same-real-directory checks."""
    raw=Path(raw_directory)
    expected=Path(expected_evaluation)
    if (raw.name!='scf_000000' or raw.resolve(strict=True)!=(expected/'scf_000000').resolve(strict=True)
            or not raw.resolve(strict=True).is_relative_to(Path(allowed_root).resolve(strict=True))):
        raise ValueError('recorded raw directory is not the exact expected native calculation')
    return raw.parent


def main():
    if OUT.exists():raise FileExistsError('fresh external audit namespace required')
    bundle=ROOT/'bundle'
    if sha(bundle/'bundle_files.json')!='d4f7a2fbf1950001b725eed0ff1d286e98825c23d62d9a070ac5bf6d4aaa00d2':
        raise ValueError('immutable execution bundle required')
    files=json.loads((bundle/'bundle_files.json').read_bytes())
    if len(files)!=9 or any(sha(bundle/n)!=h for n,h in files.items()):raise ValueError('bundle drift')
    spec=importlib.util.spec_from_file_location('original_matched_execution_worker',bundle/'worker.py')
    w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w)
    config=json.loads((bundle/'execution.json').read_bytes())
    m,c,point,old,data=w.validate(config,bundle);w.verify_runtime(m)
    accepted=json.loads((ROOT/'accepted.json').read_bytes())
    if accepted['job_id']!='28866132':raise ValueError('exact accepted matched-arm handle required')
    accounting=subprocess.run(['sacct','-n','-P','-X','-j','28862441,28865656,28866132',
        '--format=JobIDRaw,State,ExitCode,ElapsedRaw,AllocCPUS,NodeList,Partition'],
        check=True,capture_output=True,text=True).stdout.strip()
    queue=subprocess.run(['squeue','-h','-j','28862441,28865656,28866132',
        '-o','%i|%100j|%T|%P|%C|%D|%M|%R'],check=True,capture_output=True,text=True).stdout.strip()
    directory=ROOT/'arm/new_evaluations'
    paths=sorted(directory.glob('eval-*/result.json'))
    if not paths:raise RuntimeError('no completed new SCF yet; preserve live job, no restart')
    path=paths[0];row=json.loads(path.read_bytes());x=np.array(row['coordinates'])
    digest=hashlib.sha256(x.tobytes()).hexdigest()
    if (path.parent.name!=f'eval-000000-{digest[:12]}' or row['coordinate_sha256']!=digest
            or row['contract_sha256']!=old[0]['contract_sha256'] or row['status']!='complete'):
        raise ValueError('first complete record or physical contract differs')
    atoms=data['chart'].to_atoms(x)
    atoms.calc=SinglePointCalculator(atoms,**m.faces.study.audited_results(path.parent/'scf_000000'))
    recorded_directory=recorded_evaluation_directory(row['validation']['raw_directory'],path.parent,ROOT/'arm')
    native=m.faces.study.result_validator(0,atoms,recorded_directory)
    e,g=data['chart'].energy_and_gradient(x,atoms)
    internal=data['basis'].T@x
    checks=dict(native_validation_equal=native==row['validation'],energy_equal=e==row['enthalpy_eV'],
        gradient_equal=bool(np.allclose(g,row['gradient'],atol=2e-12,rtol=0)),
        stress_equal=bool(np.array_equal(atoms.get_stress(voigt=False),row['stress_eV_per_A3'])),
        Q_equal=bool(np.allclose(data['axes'].T@internal,[-.2,-.1],atol=2e-12,rtol=0)),
        prescribed_control_exact=x[-1]==0,
        roundtrip_original_guard=m.faces.canary.geometry_guard(data['chart'],data['basis'],internal))
    checks={k:bool(v) for k,v in checks.items()}
    if not all(checks.values()):
        print(json.dumps(dict(checks=checks,internal_roundtrip_norm_A=float(np.linalg.norm(internal)),
            physical_full_norm_A=float(np.linalg.norm(x)),maximum_gradient_difference=float(np.max(np.abs(g-row['gradient']))),
            Q=(data['axes'].T@internal).tolist(),native=native,recorded_validation=row['validation'])),flush=True)
        raise ValueError('native EFS/full gradient/fixed Q/control/original geometry differs')
    probe=RUNS/'20261015-varneb-G3-feasible-preflight-E170-r1/native_matched_probe.json'
    if sha(probe)!='2c22c66c58511c51d768271c17f0b6cadec7dfb6781af5106cfbd0822d5ce7c7':
        raise ValueError('original zero-DFT novel-proposal proof required')
    expected=json.loads(probe.read_bytes())['first_novel_full_A']
    gi=data['basis'].T@g;gr=gi-data['axes']@(data['axes'].T@gi)
    result=dict(status='first_new_native_SCF_audited_not_optimizer_outcome',
        observed_CST=datetime.now(timezone(timedelta(hours=8))).isoformat(),job_id='28866132',
        sacct=accounting,queue=queue,first_result_sha256=sha(path),validation=native,
        first_new_coordinate_exactly_matches_E170_novel_proposal=np.array_equal(x,expected),
        recorded_raw_path=row['validation']['raw_directory'],canonical_raw_path=str((path.parent/'scf_000000').resolve()),
        raw_directory_identity_checked_before_exact_validation_comparison=True,
        full_physical_gradient_eV_A=g.tolist(),released_gradient_norm_eV_A=float(np.linalg.norm(gr)),
        first_new_energy_eV_cell=e,total_internal_radius_A=float(np.linalg.norm(internal)),
        completed_new_SCFs_at_observation=len(paths),reserved_attempts_at_observation=len(list(directory.glob('attempt-*.json'))),
        actual_backend_failure_paths=[str(p) for p in directory.glob('eval-*/failure.json')],
        final_arm_summary_present=(ROOT/'arm/arm.json').is_file(),original_runtime_files_rechecked=380,
        additional_DFT_or_scheduler_writes=0,original_source_inputs_changed=False,
        physical_reserved_geometry_generated_or_read=False,
        convergence_stability_or_speedup_certified=False)
    OUT.mkdir()
    with (OUT/'first_native.json').open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(result))


if __name__=='__main__':main()
