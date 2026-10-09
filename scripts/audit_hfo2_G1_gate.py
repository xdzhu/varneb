"""Review the registered G1 candidate network before the single G2 canary.

This is a candidate-connectivity gate, not a TS/global-stability certificate.
All observations and sampling checks are replayed; no DFT is launched.
"""
from __future__ import annotations

import argparse
import importlib
import json
from pathlib import Path
import tempfile

import numpy as np

from examples.hfo2_fixed_input_factory import CONTRACT
from scripts.analyze_hfo2_channel_network import read_evaluated_observation
from scripts.analyze_hfo2_G1_peak_sampling import analyze as five_samples
from scripts.audit_hfo2_static_replica import audited_results, sha256
from scripts.hfo2_switching_path_polarization import save
from scripts.plot_hfo2_G1_terminal_network import build_data as terminal_data
from scripts.plot_hfo2_sampling_bridge import build_data as preserving_samples
from scripts.prepare_hfo2_reversing_peak_sampling import analyze as reversing_samples

CASE = Path("benchmarks/hfo2_channels/20261008")
EXPECTED = (
    ("PO_to_T", "b0d42fa6797986b038366f6856ba0943d2e2f4930356fd958ac638d1a2bef844", "P4_2/nmc"),
    ("PO_to_M", "d2ad5af18b18820cc92a3176d727feda13881d80dd0c91d682ee10805566766f", "P2_1/c"),
    ("PO_flip_T_pattern_preserving", "2058e4045e57ddfa23c95d551593aa861e6ccca169b9e98d3f12bf112bd5b229", "Pca2_1"),
    ("PO_flip_T_pattern_reversing", "23c1f7e3799ab0f1b388bd2791418c400f26eccf5679a3353c05dff326ef3f08", "Pca2_1"),
)


def validate_topology(channels):
    if len(channels) != 4:
        raise ValueError("four registered candidate edges required")
    for channel, (name, digest, final_group) in zip(channels, EXPECTED):
        if (channel["name"] != name or channel["source_observation_sha256"] != digest
                or not channel["ordinary_residual_passed"]
                or channel["source_NEB_fmax_eV_A"] > .10):
            raise ValueError("registered terminal ordinary candidate changed")
        rows = channel["rows"]
        if (rows[0]["energy_relative_PO_meV_fu"] != 0.
                or rows[0]["s_normalized"] != 0. or rows[-1]["s_normalized"] != 1.):
            raise ValueError("common initial well and ordered complete arc required")
        for row in rows:
            values = [row[key] for key in ("Q_pattern_x_A", "Q_pattern_y_A", "Q_pattern_z_A",
                "triplet_captured_squared_norm_fraction", "T_ordered_chart_translation_free_atomic_RMS_A",
                "T_Green_strain_Frobenius_norm")]
            if not np.isfinite(values).all():
                raise ValueError("finite mode/strain/reconstruction diagnostics required")
        for row, expected_group in ((rows[0], "Pca2_1"), (rows[-1], final_group)):
            if [s["symbol"] for s in row["structure_audit"]["symmetry_sweep"]] != [expected_group]*3:
                raise ValueError("endpoint identity not consistent at registered tolerances")


def audit(repository, *, raw_on_hf=False):
    repository = Path(repository)
    root = repository / CASE
    data = terminal_data(repository)
    validate_topology(data["channels"])
    spec_path = root / "reversing_peak_sampling_20261009/network_specification.json"
    spec = json.loads(spec_path.read_text())
    endpoints, barriers, fresh_raw = [], [], 0
    for entry, channel in zip(spec["channels"], data["channels"]):
        images, observation, digest = read_evaluated_observation(spec_path.parent / entry["observation"])
        if digest != channel["source_observation_sha256"]:
            raise ValueError("terminal and numeric replay source differ")
        view = images[::-1] if entry["reverse"] else images
        energies = np.array([a.get_potential_energy() for a in view])
        forward, reverse = (energies.max()-energies[[0, -1]])*250
        difference = (energies[-1]-energies[0])*250
        if abs(forward-reverse-difference) > 1e-8:
            raise ValueError("forward/reverse/common-well identity failed")
        barriers.append({"name": channel["name"], "terminal_sampled_forward_meV_fu": float(forward),
                         "terminal_sampled_reverse_meV_fu": float(reverse),
                         "endpoint_difference_meV_fu": float(difference)})
        for label, atoms in (("initial", view[0]), ("final", view[-1])):
            force = float(np.linalg.norm(atoms.get_forces(), axis=1).max())
            stress = float(np.abs(atoms.get_stress()).max()*1602.176634)
            if force >= .03 or stress >= 2.:
                raise ValueError("free-cell endpoint physical screen failed")
            endpoints.append({"channel": channel["name"], "end": label,
                              "max_atomic_force_eV_A": force, "max_stress_component_kbar": stress})
        if raw_on_hf:
            for atoms, point in zip(images, observation["raw_image_evaluations"]):
                directory = Path(point["raw_source"])
                if ({n: sha256(directory/n) for n in CONTRACT} != CONTRACT
                        or sha256(directory/"OUT.ABACUS/running_scf.log") != point["raw_log_sha256"]):
                    raise ValueError("original raw physical bytes/log changed")
                actual = audited_results(directory)
                expected = {"energy": atoms.get_potential_energy(), "forces": atoms.get_forces(), "stress": atoms.get_stress()}
                if any(not np.allclose(actual[k], expected[k], atol=1e-12, rtol=0) for k in expected):
                    raise ValueError("original raw E/F/stress changed")
                fresh_raw += 1
    with tempfile.TemporaryDirectory(prefix="hfo2-G1-replay-") as folder:
        t_m = five_samples(repository, root/"G1_peak_sampling_20261009/completed_HF", Path(folder)/"five.json")
        reverse = reversing_samples(repository, root/"reversing_peak_sampling_20261009/completed_HF/prepared", Path(folder)/"reverse.json")
    preserve = preserving_samples(root/"switching_converged_update_20261009_1255/preserving_step69",
                                   root/"preserving_sampling_bridge_20261009/completed_HF")
    if not (t_m["channels"]["T_to_PO"]["ordinary_residual_passed"]
            and preserve["reconstructed_band"]["ordinary_residual_passed"] and reverse["ordinary_residual_passed"]):
        raise ValueError("sampling-reconstructed candidate failed; review before canary")
    checker = importlib.import_module("benchmarks.hfo2_channels.20261008.switching_path_polarization_20261009.check_export")
    berry = checker.audit(repository, root/"switching_path_polarization_20261009/completed_HF")
    if not all(c["sampled_branch"]["status"] == "conditional_sampled_lift"
               and all(link["unique"] for link in c["sampled_branch"]["links"]) for c in berry["channels"]):
        raise ValueError("longitudinal sampled branch unresolved")
    gamma_path = root/"gamma_analysis/audit.json"
    gamma = json.loads(gamma_path.read_text())["families"]
    optical = {k: {"negative_optical_indices": gamma[k]["negative_optical_indices"],
                   "minimum_positive_optical_frequency_THz": min(gamma[k]["frequencies_thz"][3:])}
               for k in ("T_d0.01", "T_d0.02")}
    return {"status": "registered_G1_candidate_review_passed_canary_only",
            "candidate_G1_gate_passed": True, "new_DFT_calls": 0, "raw_on_HF": raw_on_hf,
            "fresh_raw_terminal_records": fresh_raw, "existing_terminal_records": 40,
            "physical_contract": CONTRACT, "common_PO_energy_eV_cell": -9783.249675811956,
            "endpoints": endpoints, "terminal_barriers": barriers,
            "sampling": {"peak_SCF_calls": 13, "cap": 14,
                "T_union_peak_meV_fu": t_m["channels"]["T_to_PO"]["highest_union_sample_meV_fu"],
                "M_refined_12_image_peak_meV_fu": barriers[1]["terminal_sampled_forward_meV_fu"],
                "preserving_union_peak_meV_fu": max(preserve["old_max_meV_fu"], preserve["new_max_meV_fu"]),
                "reversing_union_peak_meV_fu": reverse["highest_union_sample_meV_fu"]},
            "T_existing_Gamma_atomic_evidence": optical, "T_Gamma_audit_sha256": sha256(gamma_path),
            "R3_conditional_lifts": berry["channels"], "analysis_driver_sha256": sha256(Path(__file__)),
            "next_allowed_experiment": "one zero-strain PO_plus four-step BFGS clamped canary, max5 fresh SCFs",
            "G2_matrix_or_holdout_automatically_submitted": False,
            "TS_global_escape_stability_spontaneous_P_or_distinct_winding_certified": False,
            "limitations": ["G1 candidate connectivity and source/endpoint/sampling/property review, not G3 stationary bottleneck certification",
                "T Gamma optical evidence checks only its fixed-cell atomic Gamma subspace, not all q or joint atomic/cell stability",
                "Pbcn central branch stability unmeasured; M/T are not exhaustive escape coverage",
                "sampling union maxima are not error-bounded continuous activation barriers",
                "two geometric switching candidates are not certified distinct MEP or winding sectors",
                "new mechanical conditions require new unforced endpoint relaxation and branch/phase audit"]}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repository", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--raw-on-HF", action="store_true")
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError("refusing previous G1 review")
    result = audit(a.repository, raw_on_hf=a.raw_on_HF)
    save(a.output, result)
    print(json.dumps({k: result[k] for k in ("status", "new_DFT_calls", "fresh_raw_terminal_records")}))
