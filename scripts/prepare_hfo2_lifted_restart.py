"""Explicitly re-lift a complete stopped pilot snapshot, reusing exact SCFs.

No DFT, implicit atom mapping, physical input change or Slurm submission.
The original run, partial calls and source archive remain unchanged. A fresh
optimizer is required because the old tangents/springs used a broken lift.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from ase.io import read, write

from examples.hfo2_fixed_input_factory import CONTRACT, same_ordered_geometry
from scripts.audit_hfo2_static_replica import audited_results, sha256
from scripts.prepare_hfo2_reference_variants import write_clean_poscar
from vcneb import minimum_image_path_lift, validate_periodic_path_lift, validate_path_geometry


def exact_cached_source(image, directory):
    """Find identical ordered periodic geometry, not nearest-neighbor matching."""
    cache = directory/'seed_cache_audit.json'
    candidates=[]
    if cache.exists():
        candidates.append(Path(json.loads(cache.read_text())["raw_source"]))
    candidates += sorted((p for p in directory.glob('scf_*') if (p/'call_audit.json').is_file()),reverse=True)
    for source in candidates:
        if not same_ordered_geometry(image,read(source/'STRU',format='abacus')):
            continue
        if any(sha256(source/name)!=digest for name,digest in CONTRACT.items()):
            raise ValueError('matching cached SCF violates the fixed physical contract')
        raw=audited_results(source)
        call=source/'call_audit.json'
        if call.exists():
            recorded=json.loads(call.read_text())
            if (any(sha256(source/n)!=h for n,h in recorded['input_sha256'].items())
                    or sha256(source/'OUT.ABACUS/running_scf.log')!=recorded['raw_log_sha256']):
                raise ValueError('matching SCF raw evidence changed')
        return source,raw
    raise ValueError(f'no complete exact ordered SCF for snapshot image in {directory}')


def prepare(workdir,snapshot,output,source_job_id):
    if output.exists():
        raise FileExistsError('refusing existing restart namespace')
    if not snapshot.resolve().is_relative_to((workdir/'snapshots').resolve()):
        raise ValueError('snapshot must belong to the declared stopped workdir')
    files=sorted(snapshot.glob('POSCAR_*'))
    if len(files) not in (9,10) or [p.name for p in files]!=[f'POSCAR_{i:02d}' for i in range(len(files))]:
        raise ValueError('complete nine/ten-image snapshot required')
    original=[read(p,format='vasp') for p in files]
    images,lift=minimum_image_path_lift(original)
    if not np.any(lift['integer_lattice_shifts_by_image_atom']):
        raise ValueError('snapshot needs no lift repair; use ordinary resume instead')
    validate_periodic_path_lift(images)
    geometry=validate_path_geometry(images,minimum_distance=1.6,maximum_deformation=.25)
    sources,audits=[],[]
    for i,image in enumerate(images):
        source,raw=exact_cached_source(image,workdir/f'image_{i:04d}')
        sources.append(str(source))
        audits.append({'image_index':i,'raw_source':str(source),
            'energy_eV_cell':float(raw['energy']),
            'raw_log_sha256':sha256(source/'OUT.ABACUS/running_scf.log'),
            'input_sha256':{n:sha256(source/n) for n in (*CONTRACT,'STRU')},
            'snapshot_POSCAR_sha256':sha256(files[i]),
            'snapshot_to_raw_ordered_periodic_match':True})
    output.mkdir(parents=True)
    write(output/'seed.traj',images)
    write_clean_poscar(output/'initial.vasp',images[0])
    write_clean_poscar(output/'final.vasp',images[-1])
    (output/'factory_parameters.json').write_text(json.dumps({
        'source_directory':sources[0],'seed_static_directories':sources},indent=2)+'\n')
    report={'purpose':'explicit_periodic_lift_repair_of_stopped_health_pilot',
        'source_job_id':source_job_id,'source_workdir':str(workdir),'complete_snapshot':str(snapshot),
        'n_total_images':len(images),'n_fixed_endpoints':2,'n_active_images':len(images)-2,
        'fmax_eV_A':.10,'climb':False,'pressure_GPa':0,'optimizer':'FIRE',
        'maxstep_A':.02,'spring_eV_A2':.2,'first_optimization_step_cap':10,
        'physical_inputs_changed':False,'periodic_lift_audit':lift,'initial_geometry_gate':geometry,
        'seed_cache_audits':audits,'old_optimizer_state_reused':False,
        'DFT_executed':False,'job_submitted':False,
        'limitations':'explicit short-adjacent-step convention; old forces are not convergence or acceleration evidence'}
    (output/'manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('workdir','snapshot','output'):
        parser.add_argument(f'--{key}',type=Path,required=True)
    parser.add_argument('--source-job-id',required=True)
    args=parser.parse_args()
    report=prepare(args.workdir,args.snapshot,args.output,args.source_job_id)
    print(json.dumps({'prepared':True,'n_total_images':report['n_total_images'],
                      'all_seed_SCFs_reused':True,'DFT_executed':False,'job_submitted':False}))


if __name__=='__main__':
    main()
