"""Make a separate audited snapshot from exact-state cached evaluations.

This is an IO-only archival operation, not path optimization, SCF or a
convergence certificate. Never writes into or changes the production inputs.
Requires the original fully evaluated chain as independent corroboration.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read, write

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from vcneb.executor import ImageEvaluation, ThreadedCalculatorExecutor
from vcneb.vasp_contract import parameter_digest


def file_sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def physics_namespace(preflight):
    command = preflight["command_line"]
    flag = "--image-cache-namespace"
    user = command[command.index(flag)+1] if flag in command else None
    return parameter_digest({"parameters": preflight["calculator_parameters"],
        "inputs": {k: v["sha256"] for k,v in preflight["licensed_input_fingerprints"].items()},
        "vca": preflight.get("vca"), "user_namespace": user, "input_contract_version": 1})


def same_geometry(a, b):
    return all(np.array_equal(x,y) for x,y in ((a.positions,b.positions), (a.cell.array,b.cell.array),
                                               (a.numbers,b.numbers), (a.pbc,b.pbc)))


def restore(snapshot, reference, preflight_path, reference_preflight_path, cache_dir, output):
    snapshot, reference, preflight_path, reference_preflight_path, cache_dir, output = map(Path,
        (snapshot, reference, preflight_path, reference_preflight_path, cache_dir, output))
    if output.resolve() in (snapshot.resolve(), reference.resolve()) or output.exists() or output.with_suffix(".json").exists():
        raise FileExistsError("independent new output required; source/reference cannot be replaced")
    hashes = {str(p):file_sha(p) for p in (snapshot,reference,preflight_path,reference_preflight_path)}
    preflight = json.loads(preflight_path.read_text())
    ref_preflight = json.loads(reference_preflight_path.read_text())
    namespace = physics_namespace(preflight)
    if namespace != physics_namespace(ref_preflight):
        raise ValueError("reference physics namespace differs")
    meta_path = cache_dir / "cache_metadata.json"
    meta = json.loads(meta_path.read_text())
    if meta != {"format_version":1,"namespace":namespace}:
        raise ValueError("cache physics namespace differs")
    cache_meta_sha = file_sha(meta_path)
    images, corroboration = read(snapshot,index=":"), read(reference,index=":")
    if len(images) != int(preflight["n_images"]) or len(corroboration) != len(images):
        raise ValueError("a single complete chain is required")
    records, restored = [], []
    for index, (image, ref) in enumerate(zip(images,corroboration)):
        if not same_geometry(image,ref):
            raise ValueError(f"reference geometry is not byte-exact for image {index}")
        if ref.calc is None or not all(k in ref.calc.results for k in ("energy","forces","stress")):
            raise ValueError(f"reference image {index} lacks complete evaluated results")
        verified = ImageEvaluation(ref.get_potential_energy(), ref.get_forces(), ref.get_stress(voigt=False)).validate(image_index=index,n_atoms=len(ref))
        existing = image.calc is not None and all(k in image.calc.results for k in ("energy","forces","stress"))
        if index in (0,len(images)-1):
            if not existing:
                raise ValueError("endpoint must already contain its original static result")
            value = ImageEvaluation(image.get_potential_energy(),image.get_forces(),image.get_stress(voigt=False)).validate(image_index=index,n_atoms=len(image))
            record = {"image_index":index,"source":"existing_endpoint_static"}
        else:
            digest = ThreadedCalculatorExecutor._image_cache_digest(index,image)
            path = cache_dir / f"image_{index:04d}_{digest}.npz"
            if not path.is_file():
                raise FileNotFoundError(f"no exact-state cache for image {index}: {path}")
            cache_sha = file_sha(path)
            with np.load(path,allow_pickle=False) as data:
                value = ImageEvaluation(float(np.asarray(data["energy"]).reshape(())),data["forces"],data["stress"]).validate(image_index=index,n_atoms=len(image))
            if file_sha(path) != cache_sha:
                raise RuntimeError("cache file changed during read")
            if existing and (image.get_potential_energy()!=value.energy or not np.array_equal(image.get_forces(),value.forces)
                             or not np.array_equal(image.get_stress(voigt=False),value.stress)):
                raise ValueError("existing snapshot results conflict with exact-state cache")
            record = {"image_index":index,"source":"exact_state_cache","cache_key":digest,"cache_path":str(path),"cache_sha256":cache_sha,"filled_missing_results":not existing}
        if value.energy != verified.energy or not np.array_equal(value.forces,verified.forces) or not np.array_equal(value.stress,verified.stress):
            raise ValueError(f"cache/existing results disagree with independent reference for image {index}")
        copy = image.copy()
        copy.calc = SinglePointCalculator(copy,energy=value.energy,forces=value.forces,stress=value.stress)
        restored.append(copy)
        records.append(record)
    for path, digest in hashes.items():
        if file_sha(Path(path)) != digest:
            raise RuntimeError("source changed during archival check")
    if file_sha(meta_path) != cache_meta_sha:
        raise RuntimeError("cache metadata changed during archival check")
    output.parent.mkdir(parents=True,exist_ok=True)
    # Reserve exclusively, after validating the entire chain. Never replace an
    # existing artifact. Output is unrelated to the active production files.
    with output.open("xb") as handle:
        write(handle,restored,format="traj")
    saved = read(output,index=":")
    if not all(same_geometry(a,b) for a,b in zip(images,saved)):
        raise RuntimeError("archival roundtrip changed geometry")
    report = {"status":"exact_state_archival_copy_verified","not_a_convergence_certificate":True,
              "DFT_calls":0,"geometry_modified":False,"source_files_sha256":hashes,
              "physics_namespace":namespace,"cache_metadata_sha256":cache_meta_sha,
              "output":str(output),"output_sha256":file_sha(output),"n_images":len(saved),"images":records}
    with output.with_suffix(".json").open("x") as handle:
        json.dump(report,handle,indent=2)
        handle.write("\n")
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ("snapshot","reference","preflight","reference-preflight","cache-dir","output"):
        p.add_argument("--"+arg,required=True,type=Path)
    a=p.parse_args()
    r=restore(a.snapshot,a.reference,a.preflight,a.reference_preflight,a.cache_dir,a.output)
    print(json.dumps({k:r[k] for k in ("status","DFT_calls","geometry_modified","n_images","output_sha256")},indent=2))


if __name__=="__main__":
    main()
