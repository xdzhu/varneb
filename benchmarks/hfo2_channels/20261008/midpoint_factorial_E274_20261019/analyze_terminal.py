"""Native four-state energy contrasts, not a relaxed-path or forecast update."""
from __future__ import annotations
import argparse
import csv
import importlib.util
import json
from pathlib import Path
import sys
import tarfile
import tempfile

import numpy as np

ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT))
from examples.hfo2_fixed_input_factory import CONTRACT,read_fixed_hfo2_stru,same_ordered_geometry
from scripts.audit_hfo2_static_replica import audited_results,sha256

CASE=Path(__file__).resolve().parent
PARENT=CASE.with_name('midpoint_reverse_E269_20261018')
PIN='d9cc4c2b81ecf7dabb542c943913a75783e5e0a1d472c0a47c591588bbb4f6ca'
JOB='28950839'


def load(name):
    spec=importlib.util.spec_from_file_location('E274_original_'+name,PARENT/f'{name}.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def restore(destination):
    receipt=json.loads((CASE/'native_terminal_receipt.json').read_bytes());archive=CASE/'native_terminal.tar.gz'
    if (sha256(archive)!=receipt['archive_sha256'] or receipt['job_id']!=JOB or receipt['bundle_pin']!=PIN
            or receipt['new_SCF_calls']!=4 or receipt['new_DFT_calls_for_analysis']!=0
            or receipt['scheduler_mutations']!=0):
        raise ValueError('complete original four-call native terminal required')
    with tarfile.open(archive) as tar:
        for member in tar.getmembers():
            p=Path(member.name)
            if (p.is_absolute() or '..' in p.parts or not p.parts or p.parts[0]!='native_terminal'
                    or not(member.isdir() or member.isfile())
                    or not(destination/p).resolve().is_relative_to(destination.resolve())):
                raise ValueError('unsafe native terminal member')
        tar.extractall(destination)
    folder=destination/'native_terminal';files=json.loads((folder/'files.json').read_bytes())
    if len(files)!=receipt['exported_files']:raise ValueError('complete portable file count differs')
    for name,digest in files.items():
        p=folder/name
        if not p.resolve().is_relative_to(folder.resolve()) or sha256(p)!=digest:
            raise ValueError('native terminal bytes changed')
    return folder,receipt


def analyze():
    cohort=load('analyze_completed');preparer=load('prepare_factorial')
    with tempfile.TemporaryDirectory(prefix='varneb_E274_native_') as tmp:
        temporary=Path(tmp);folder,receipt=restore(temporary)
        complete=json.loads((folder/'complete.json').read_bytes())
        prep=json.loads((folder/'preparation.json').read_bytes())
        proof=json.loads((folder/'preflight.json').read_bytes())
        row=json.loads((folder/'scheduler.json').read_bytes())['row']
        if (row[0]!=JOB or row[1:7]!=['COMPLETED','0:0','32','hfacnormal01','iai806','1']
                or complete['job_id']!=JOB or complete['bundle_pin']!=PIN or complete['new_SCF_calls']!=4
                or len(complete['records'])!=4
                or prep!=json.loads((CASE/'prepared/preparation.json').read_bytes())
                or sha256(folder/'preflight.json')!=receipt['native_preflight_sha256']
                or sha256(CASE/'native_preflight.json')!=receipt['native_preflight_sha256']
                or proof['bundle_pin']!=PIN or proof['new_DFT_calls']!=0):
            raise ValueError('original completed job, preparation and native input proof must agree')
        native={};seconds=[]
        for record,deck in zip(complete['records'],prep['decks']):
            channel,corner=record['channel'],record['corner'];raw=folder/channel/corner
            if (channel,corner)!=(deck['channel'],deck['corner']) or (channel,corner) in native:
                raise ValueError('all four distinct ordered crossover records required')
            call=json.loads((raw/'call_audit.json').read_bytes());parsed=audited_results(raw)
            native_proof=next(r for r in proof['decks'] if (r['channel'],r['corner'])==(channel,corner))
            if (any(call['input_sha256'][n]!=h for n,h in CONTRACT.items())
                    or any(sha256(raw/n)!=call['input_sha256'][n] for n in ('INPUT','KPT','STRU'))
                    or call['input_sha256']!=record['input_sha256']
                    or call['input_sha256']['STRU']!=native_proof['input_sha256']['STRU']
                    or sha256(raw/'OUT.ABACUS/running_scf.log')!=record['raw_log_sha256']
                    or call['raw_log_sha256']!=record['raw_log_sha256']
                    or any(not np.allclose(parsed[k],record['results'][k],atol=1e-12,rtol=0)
                           or not np.allclose(parsed[k],call['results'][k],atol=1e-12,rtol=0) for k in parsed)
                    or not same_ordered_geometry(read_fixed_hfo2_stru(raw/'STRU'),
                        read_fixed_hfo2_stru(CASE/'prepared'/deck['path']/'STRU'))):
                raise ValueError('native frozen geometry/input/energy/forces/stress differs')
            if not np.isfinite(record['seconds']) or record['seconds']<=0 or record['seconds']!=call['elapsed_seconds']:
                raise ValueError('positive actual timing required')
            native[channel,corner]=parsed;seconds.append(record['seconds'])
        if len(native)!=4 or len(list(folder.glob('*/*/call_audit.json')))!=4:
            raise ValueError('exactly four fresh native SCFs required')
        exported,_=cohort.restore(temporary);rows=[];corners=[]
        for channel in ('preserving','T'):
            source=next(r for r in prep['channels'] if r['channel']==channel)
            energies={}
            for label,path in (('A',ROOT/preparer.OLD_RESULTS/channel),('B',exported/channel)):
                atoms,raw_record,raw,provenance=preparer.source_peak(path)
                if provenance!=source['sources'][label]:raise ValueError('original diagonal source changed')
                if (any(sha256(raw/n)!=raw_record['input_sha256'][n] for n in ('INPUT','KPT','STRU'))
                        or sha256(raw/'OUT.ABACUS/running_scf.log')!=raw_record['raw_log_sha256']
                        or not same_ordered_geometry(atoms,read_fixed_hfo2_stru(raw/'STRU'))):
                    raise ValueError('native diagonal geometry, physical input or raw log differs')
                parsed=audited_results(raw)
                if parsed['energy']!=raw_record['energy_eV_cell']:raise ValueError('native diagonal energy differs')
                energies[label*2]=parsed['energy']
            energies.update({c:native[channel,c]['energy'] for c in ('AB','BA')})
            contrast=preparer.contrasts(energies)
            if not np.isclose(contrast['diagonal_change_meV_fu'],contrast['internal_at_A_cell_meV_fu']
                    +contrast['cell_at_A_internal_meV_fu']+contrast['finite_interaction_meV_fu'],atol=1e-8,rtol=0):
                raise ValueError('finite interaction closure differs')
            rows.append(dict(channel=channel,**contrast))
            corners.extend(dict(channel=channel,corner=c,energy_eV_cell=energies[c],
                relative_to_AA_meV_fu=(energies[c]-energies['AA'])*250,new_SCF=c in ('AB','BA'))
                for c in ('AA','AB','BA','BB'))
        proxy=sum(seconds)*32/3600;allocated=float(row[7])*32/3600
        if (not np.isclose(proxy,receipt['native_proxy_core_hours'],atol=1e-10,rtol=0)
                or not np.isclose(allocated,receipt['allocated_core_hours'],atol=1e-10,rtol=0)):
            raise ValueError('all actual native and allocation costs must agree')
    return dict(status='four_measured_frozen_crossovers_not_relaxed_barriers',rows=rows,corners=corners,
        source_job=JOB,bundle_pin=PIN,native_archive_sha256=receipt['archive_sha256'],
        all_four_new_native_EFS_reparsed_locally=True,native_diagonal_energies_replayed_locally=True,
        full_pseudo_orbital_bytes_rechecked_locally=False,new_SCF_calls=4,
        native_proxy_core_hours=proxy,allocated_core_hours=allocated,earlier_costs_additional=True,
        new_DFT_calls_for_analysis=0,old_labels_and_forecasts_unchanged=True,
        unique_elastic_or_phonon_energy_partition=False,new_relaxed_path_or_TS_certified=False)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args()
    result=analyze();a.output.mkdir(parents=True,exist_ok=False)
    (a.output/'analysis.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    for name,rows in (('contrasts',result['rows']),('corner_energies',result['corners'])):
        with (a.output/f'{name}.csv').open('w',encoding='utf-8',newline='') as handle:
            writer=csv.DictWriter(handle,list(rows[0]));writer.writeheader();writer.writerows(rows)
    print(json.dumps(result,allow_nan=False))
