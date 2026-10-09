"""Prepare one audited PO+ geometry continuation, zero new DFT.

No extra phase/strain/channel is generated. A fresh BFGS Hessian is used;
the identical terminal clamped SCF is reused, not a free-cell evaluation.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np
from ase.io import read

from examples.hfo2_fixed_input_factory import read_fixed_hfo2_stru, same_ordered_geometry
from scripts.audit_hfo2_clamped_canary import audit
from scripts.audit_hfo2_static_replica import sha256
from scripts.hfo2_switching_path_polarization import save
from scripts.prepare_hfo2_reference_variants import write_clean_poscar
from scripts.relax_clamped_ase_endpoint import load_seed
from vcneb import endpoint_structure_record

CASE = Path("benchmarks/hfo2_channels/20261008")
CANARY_HF = "/public/home/iai806/abacus/agent-runs/20261009-varneb-clamped-canary-E045-r1"
ORIGINAL_CANARY_AUDIT_SHA256 = "f2e0721caddde5c868235f776eabddc3acca76e690f3095c85840f1621c09ce5"


def prepare(repository, canary_endpoint, output, *, full_physical_bytes=True):
    if output.exists():
        raise FileExistsError("refusing an existing continuation namespace")
    case = repository/CASE
    seed_manifest = case/"clamped_endpoint_seeds/strain_0000/PO_plus/endpoint_seed.json"
    _, boundary, original = load_seed(seed_manifest)
    parent_audit = canary_endpoint/"canary_audit.json"
    if sha256(parent_audit) != ORIGINAL_CANARY_AUDIT_SHA256:
        raise ValueError("registered original canary audit changed")
    review = audit(canary_endpoint, seed_manifest, full_physical_bytes=full_physical_bytes)
    if (review["new_SCF_calls"] != 5 or review["endpoint_status"] != "step_limit"
            or review["endpoint_converged"]):
        raise ValueError("only the actual four-step, still-unconverged canary is continued")
    if any([s["symbol"] for s in row["structure_audit"]["symmetry_sweep"]] != ["Pca2_1"]*3
           for row in review["rows"]):
        raise ValueError("canary PO identity changed; review rather than force restoring it")
    last = review["rows"][-1]
    static = canary_endpoint/"calculator/image_0000/scf_000004"
    atoms = read_fixed_hfo2_stru(static/"STRU")
    boundary.validate_images([atoms])
    output.mkdir(parents=True, exist_ok=False)
    path = output/"POSCAR.seed"
    write_clean_poscar(path, atoms)
    restored = read(path, format="vasp")
    if not same_ordered_geometry(atoms, restored):
        raise ValueError("continuation geometry changed in POSCAR round trip")
    boundary.validate_images([restored])
    distances = restored.get_all_distances(mic=True)
    np.fill_diagonal(distances, np.inf)
    manifest = {**original, "seed_file": path.name, "seed_sha256": sha256(path),
        "source_sha256": last["input_sha256"]["STRU"],
        "ordered_geometry": endpoint_structure_record(restored),
        "minimum_distance_A": float(distances.min()), "volume_A3": float(restored.get_volume()),
        "continuation": {"parent_job": "28446324", "parent_BFGS_steps": 4,
            "parent_audit_sha256": sha256(parent_audit), "parent_log_sha256": last["raw_log_sha256"],
            "parent_summary_sha256": review["summary_sha256"],
            "status_meaning": "partial clamped relaxation, not a converged endpoint or virgin geometry",
            "same_boundary_and_physical_contract": True, "restore_BFGS_Hessian": False},
        "new_DFT_calls": 0}
    parameters = json.loads((case/"M_endpoint_factory_parameters.json").read_text())
    parameters.update(seed_static_directory=CANARY_HF+"/endpoint/calculator/image_0000/scf_000004",
                      seed_input_sha256=last["input_sha256"], seed_raw_log_sha256=last["raw_log_sha256"])
    save(output/"endpoint_seed.json", manifest)
    save(output/"factory_parameters.json", parameters)
    receipt = {"status": "audited_single_PO_geometry_continuation_prepared_no_DFT",
        "parent_review_full_physical_bytes_checked_here": full_physical_bytes,
        "new_DFT_calls": 0, "endpoint_seed_manifest_sha256": sha256(output/"endpoint_seed.json"),
        "factory_parameters_sha256": sha256(output/"factory_parameters.json"),
        "preparation_driver_sha256": sha256(Path(__file__)), "planned_max_BFGS_steps": 20,
        "planned_max_new_SCF_calls": 20, "identical_initial_SCF_reused": True,
        "automatic_matrix_or_holdout_submission": False}
    save(output/"preparation.json", receipt)
    return receipt


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repository", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--canary-endpoint", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--offline-omitted-basis", action="store_true")
    a = p.parse_args()
    print(json.dumps(prepare(a.repository, a.canary_endpoint, a.output,
                            full_physical_bytes=not a.offline_omitted_basis)))
