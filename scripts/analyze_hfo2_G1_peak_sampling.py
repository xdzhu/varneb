"""Replay five actual peak SCFs offline, with no new DFT or optimization."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

import numpy as np
from ase.io import read
from ase.calculators.singlepoint import SinglePointCalculator

from scripts.prepare_hfo2_G1_peak_sampling import BASE, SOURCES, PO_ENERGY, load_manifest, save
from scripts.analyze_hfo2_channel_network import read_evaluated_observation
from scripts.audit_hfo2_static_replica import audited_results, sha256
from scripts.plot_hfo2_sampling_bridge import fixed_writer_geometry_matches
from examples.hfo2_fixed_input_factory import CONTRACT
from vcneb import VCNEB


def analyze(repository, completed, output):
    if output.exists():
        raise FileExistsError("refusing existing five-point analysis")
    m = load_manifest(completed)
    summary = json.loads((completed / "summary.json").read_text())
    if (summary["status"] != "five_peak_statics_audited_not_relaxed_MEP" or summary["new_SCF_calls"] != 5
            or len(summary["points"]) != 5 or summary["automatic_restart_or_G2_submission"]):
        raise ValueError("complete five-point summary required")
    probe = (completed / "mpi_affinity.txt").read_text().splitlines()
    matches = [re.fullmatch(r"rank=(\d+) host=\S+", line) for line in probe]
    if any(match is None for match in matches):
        raise ValueError("malformed actual MPI probe")
    ranks = [int(match.group(1)) for match in matches]
    if sorted(ranks) != list(range(32)):
        raise ValueError("actual unique 32-rank probe required")
    bands, records, elapsed = {}, [], []
    for name, source in SOURCES.items():
        images, raw, digest = read_evaluated_observation(repository / BASE / source["folder"])
        if digest != source["sha256"] or raw["source_job_id"] != source["job"] or raw["snapshot_step"] != source["step"]:
            raise ValueError("frozen ordinary source changed")
        bands[name] = images
    inserted = {name: {} for name in bands}
    for point, report in zip(m["points"], summary["points"]):
        i = point["index"]
        rawdir = completed / "calculations" / f"{i:02d}" / "scf_000000"
        actual = audited_results(rawdir)
        call = json.loads((rawdir / "call_audit.json").read_text())
        hashes = json.loads((rawdir / "input_sha256.json").read_text())
        atoms = read(completed / point["geometry"], format="vasp")
        if (any(report[k] != point[k] for k in ("index", "channel", "segment", "fraction", "geometry_sha256"))
                or report["manifest_sha256"] != sha256(completed / "manifest.json")
                or report["physical_input_sha256"] != CONTRACT or {n: hashes[n] for n in CONTRACT} != CONTRACT
                or call["input_sha256"] != hashes
                or any(sha256(rawdir/n) != CONTRACT[n] for n in ("INPUT", "KPT"))
                or sha256(rawdir / "STRU") != report["input_STRU_sha256"] or hashes["STRU"] != report["input_STRU_sha256"]
                or sha256(rawdir / "OUT.ABACUS/running_scf.log") != report["raw_log_sha256"]
                or call["raw_log_sha256"] != report["raw_log_sha256"]
                or sha256(completed / point["geometry"]) != point["geometry_sha256"]
                or not fixed_writer_geometry_matches(rawdir / "STRU", atoms)
                or any(not np.allclose(actual[k], r["results"][k], rtol=0, atol=1e-12)
                       for k in ("energy", "forces", "stress") for r in (report, call))
                or abs(report["relative_energy_meV_fu"]-(actual["energy"]-PO_ENERGY)*250) > 1e-10
                or not np.isfinite(call["elapsed_seconds"]) or call["elapsed_seconds"] <= 0):
            raise ValueError("actual raw sampling provenance changed")
        atoms.calc = SinglePointCalculator(atoms, **actual)
        inserted[point["channel"]].setdefault(point["segment"][0], []).append((point["fraction"], atoms))
        records.append({"index": i, "channel": point["channel"], "fraction": point["fraction"],
                        "relative_energy_meV_fu": report["relative_energy_meV_fu"],
                        "raw_log_sha256": report["raw_log_sha256"]})
        elapsed.append(call["elapsed_seconds"])
    channels = {}
    for name, original in bands.items():
        added = []
        for i, atoms in enumerate(original):
            added.append(atoms)
            added.extend(a for _, a in sorted(inserted[name].get(i, []), key=lambda pair: pair[0]))
        band = VCNEB(added, pressure=0., k=.2, climb=False)
        residual = float(np.linalg.norm(band.get_forces().reshape(-1, 3), axis=1).max())
        old = max((a.get_potential_energy()-PO_ENERGY)*250 for a in original)
        highest = max(r["relative_energy_meV_fu"] for r in records if r["channel"] == name)
        expected = {"old_sampled_maximum_meV_fu": old, "highest_new_sample_meV_fu": highest,
                    "highest_union_sample_meV_fu": max(old, highest), "new_minus_old_meV_fu": highest-old}
        if any(abs(summary["channels"][name][k]-v) > 1e-10 for k, v in expected.items()):
            raise ValueError("summary maxima changed")
        channels[name] = {**expected, "inserted_total_images": len(added), "moving_images": len(added)-2,
                          "replayed_ordinary_fmax_eV_A": residual, "ordinary_residual_passed": residual <= .10,
                          "optimizer_steps": 0, "additional_SCF_calls": 0,
                          "forward_sampled_barrier_meV_fu": max(old, highest)-(original[0].get_potential_energy()-PO_ENERGY)*250,
                          "reverse_sampled_barrier_meV_fu": max(old, highest)-(original[-1].get_potential_energy()-PO_ENERGY)*250,
                          "endpoint_difference_meV_fu": (original[-1].get_potential_energy()-original[0].get_potential_energy())*250}
    result = {"status": "five_actual_G1_sampling_checks_not_continuous_MEP_or_TS", "channels": channels,
              "points": records, "new_SCF_calls": 5, "unique_MPI_ranks": len(ranks),
              "SCF_seconds": elapsed, "total_SCF_seconds": sum(elapsed), "SCF_core_hours": sum(elapsed)*32/3600,
              "manifest_sha256": sha256(completed / "manifest.json"), "summary_sha256": sha256(completed / "summary.json"),
              "analysis_driver_sha256": sha256(Path(__file__)), "new_DFT_calls_for_analysis": 0,
              "automatic_restart_or_G2_submission": False, "limitations": m["limitations"],
              "physical_bytes_checked_on_HF": "runtime rechecked all six files before SCF; offline INPUT/KPT bytes and raw records rechecked; repeated UPF/orbitals not exported"}
    save(output, result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repository", type=Path, default=Path("."))
    p.add_argument("--completed", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    print(json.dumps(analyze(a.repository, a.completed, a.output)["channels"], indent=2))


if __name__ == "__main__":
    main()
