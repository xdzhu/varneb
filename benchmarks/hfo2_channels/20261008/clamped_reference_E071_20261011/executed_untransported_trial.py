"""Apply the already registered reference spaces to a terminal clamped band.

Zero DFT: reparse portable native data, replay the actual clamped boundary,
and retain every E023 reference registration. Geometry coverage is neither
an energy partition nor a predictive-accuracy or saddle certificate.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read

from examples.hfo2_fixed_input_factory import CONTRACT, read_fixed_hfo2_stru, same_ordered_geometry
from scripts.analyze_hfo2_clamped_residual import analyze as replay_residual
from scripts.analyze_hfo2_cmma_path_reference import MASS, lowest_complete_doublet_indices, prepare_frame
from scripts.analyze_hfo2_parent_patterns import rotated_t_triplet
from scripts.audit_hfo2_clamped_terminal import validate_terminal
from scripts.audit_hfo2_static_replica import audited_results, sha256
from vcneb.continuous_projection import continuous_reference_coordinates, project_reference_basis
from vcneb.mode_subspaces import real_mode_subspace
from vcneb.phonons import identify_acoustic_modes
from vcneb.qe_modes import read_qe_gamma_modes

REGISTRATION_SHA256 = "c00b08287391de48903b9989462c77a41e123587656a75f610a6c228aa0ac53c"


def pinned_json(path, digest):
    if (not isinstance(digest, str) or len(digest) != 64
            or any(c not in "0123456789abcdef" for c in digest) or sha256(path) != digest):
        raise ValueError("registered evidence digest differs")
    return json.loads(path.read_text(encoding="utf-8"))


def terminal_images(case, audit_sha256):
    """Revalidate nine portable frames, not the undistributed full runtime.

    The pinned HF receipt attests its historical complete input/runtime audit.
    Locally we reparse only the portable nine native SCFs; no new claim about
    checking undistributed orbitals, all earlier SCFs or a live Slurm job.
    """
    receipt = pinned_json(case/"audit_receipt.json", audit_sha256)
    summary_path = case/"vcneb_summary.json"
    if sha256(summary_path) != receipt["terminal_summary_sha256"]:
        raise ValueError("terminal summary changed")
    if (receipt["status"] != "audited_ordinary_converged_not_TS_or_sampling_certificate"
            or receipt["ordinary_converged"] is not True
            or receipt["new_DFT_calls"] != 0 or receipt["scheduler_mutations"] != 0):
        raise ValueError("an actual ordinary-converged terminal receipt is required")
    validate_terminal(receipt["job_id"], receipt["scheduler_row"],
                      json.loads(summary_path.read_text(encoding="utf-8")), False)
    observation = case/"observations"/f"step_{receipt['terminal_step']:04d}"
    report = json.loads((observation/"observation.json").read_text(encoding="utf-8"))
    if report != receipt["terminal_observation"]:
        raise ValueError("terminal observation differs from the pinned HF audit")
    # This also checks all cached E/F/stress against the exact original
    # clamped residual, rather than replaying a fully relaxed-cell NEB.
    residual = replay_residual(observation)
    images = read(observation/"evaluated_chain.traj", index=":")
    for i, (image, row) in enumerate(zip(images, report["raw_image_evaluations"])):
        raw = observation/"raw"/f"image_{i:04d}"
        if (image.get_chemical_symbols() != ["Hf"]*4+["O"]*8
                or set(row["input_sha256"]) != {*CONTRACT, "STRU"}
                or any(row["input_sha256"][n] != h for n, h in CONTRACT.items())
                or any(sha256(raw/n) != row["input_sha256"][n] for n in ("INPUT", "KPT", "STRU"))
                or sha256(raw/"OUT.ABACUS/running_scf.log") != row["raw_log_sha256"]
                or sha256(observation/f"POSCAR_{i:02d}") != row["snapshot_POSCAR_sha256"]
                or not same_ordered_geometry(image, read_fixed_hfo2_stru(raw/"STRU"))
                or not same_ordered_geometry(image, read(observation/f"POSCAR_{i:02d}", format="vasp"))):
            raise ValueError("portable native input/log/ordered geometry differs")
        parsed = audited_results(raw)
        for key, expected in (("energy", row["energy_eV_cell"]), ("forces", row["forces_eV_A"]),
                              ("stress", row["stress_ASE_voigt_eV_A3"])):
            if not np.allclose(parsed[key], expected, atol=1e-12, rtol=0):
                raise ValueError("reparsed native E/F/stress differs")
        image.calc = SinglePointCalculator(image, **parsed)
    if not residual["ordinary_residual_pass"]:
        raise ValueError("replayed ordinary residual no longer passes")
    return images, report, residual, receipt


def analyze(case, source_root, registration_path, variants, gamma_path, *, terminal_audit_sha256):
    registered = pinned_json(registration_path, REGISTRATION_SHA256)
    t_path, po_path = variants/"T.vasp", variants/"PO.vasp"
    if any(sha256(p) != registered[k] for p, k in (
            (t_path, "production_T_SHA256"), (po_path, "production_PO_SHA256"),
            (gamma_path, "production_T_Gamma_SHA256"))):
        raise ValueError("original T/PO/Gamma reference changed")
    source_pins = {name: row["sha256"] for name, row in registered["author_source_audit"]["sources"].items()}
    source_pins["Tetragonal.vasp"] = registered["additional_T_anchor"]["SHA256"]
    if any(sha256(source_root/n) != h for n, h in source_pins.items()):
        raise ValueError("pinned author reference bytes differ")
    images, report, residual, receipt = terminal_images(case, terminal_audit_sha256)
    t, cmma = read(t_path, format="vasp"), read(source_root/"Cmma.vasp", format="vasp")
    chart = continuous_reference_coordinates(images, t)
    _, patterns, _ = rotated_t_triplet(t)
    span = real_mode_subspace(patterns.reshape(3, 36).T*np.repeat(np.sqrt(MASS), 3)[:, None], expected_rank=3)
    triplet = project_reference_basis(chart["displacements_A"], span.basis_columns, metric_weights=MASS)
    with np.load(gamma_path, allow_pickle=False) as gamma:
        if not np.array_equal(gamma["masses_amu"], MASS):
            raise ValueError("original Gamma masses differ")
        acoustic, _ = identify_acoustic_modes(gamma["eigenvectors_mass_weighted"], MASS)
        lowest = lowest_complete_doublet_indices(gamma["frequencies_thz"], np.setdiff1d(np.arange(36), acoustic))
        low = project_reference_basis(chart["displacements_A"], gamma["eigenvectors_mass_weighted"][:, lowest], metric_weights=MASS)
    modes = read_qe_gamma_modes(source_root/"matdyn.modes", MASS)
    frames = []
    for frame in registered["selected_frame_registrations"]:
        reference, basis, diagnostics = prepare_frame(modes, cmma, t,
            registered["author_source_audit"]["QE_to_published_POSCAR"], frame)
        common = project_reference_basis(chart["displacements_A"], basis, metric_weights=MASS)
        affine_chart = continuous_reference_coordinates(images, reference)
        delta = affine_chart["displacements_A"]-chart["displacements_A"]
        if not np.allclose(delta, delta[0], atol=1e-10, rtol=0):
            raise ValueError("reference origin offset is not constant over the ordered path")
        affine = project_reference_basis(affine_chart["displacements_A"], basis, metric_weights=MASS)
        frames.append(dict(frame_id=frame["id"], basis_audit=diagnostics,
            fixed_initial_integer_gauge=affine_chart["fixed_integer_lattice_shifts_by_atom"].tolist(),
            common_T_origin_fraction=common["captured_squared_norm_fraction"].tolist(),
            common_T_origin_residual_sqrt_amu_A=common["residual_norm"].tolist(),
            Cmma_affine_origin_residual_sqrt_amu_A=affine["residual_norm"].tolist(),
            Cmma_affine_origin_fraction_not_same_denominator=affine["captured_squared_norm_fraction"].tolist()))
    energies = np.array([row["energy_eV_cell"] for row in report["raw_image_evaluations"]])
    relative = (energies-energies[0])*250
    if not np.allclose(relative, report["relative_enthalpy_meV_fu"], atol=1e-9, rtol=0):
        raise ValueError("native energy profile differs from the audited P=0 chain")
    return dict(status="ordinary_clamped_band_descriptive_reference_analysis_not_prediction",
        script_sha256=sha256(Path(__file__)), terminal_audit_sha256=terminal_audit_sha256,
        original_registration_sha256=REGISTRATION_SHA256, author_source_sha256=source_pins,
        source_job_id=receipt["job_id"], snapshot_step=receipt["terminal_step"],
        trajectory_sha256=report["evaluated_chain_sha256"],
        portable_native_frames_reparsed=9, complete_six_physical_bytes_rechecked_locally=False,
        prior_HF_native_SCFs_audited=receipt["fresh_interior_SCFs"],
        mechanical_boundary=report["mechanical_boundary"], pressure_GPa=0,
        replayed_fmax_eV_A=residual["fmax_eV_A"], highest_image_index=int(np.argmax(energies)),
        relative_energy_meV_fu=relative.tolist(),
        fixed_T_chart_integer_gauge=chart["fixed_integer_lattice_shifts_by_atom"].tolist(),
        green_strains_from_original_T=chart["green_strains"].tolist(),
        T_origin_zero_displacement=triplet["zero_displacement"].tolist(),
        T_triplet_rank3_fraction=triplet["captured_squared_norm_fraction"].tolist(),
        T_triplet_residual_sqrt_amu_A=triplet["residual_norm"].tolist(),
        T_lowest_two_complete_doublets_rank4_indices=lowest.tolist(),
        T_lowest_two_complete_doublets_rank4_fraction=low["captured_squared_norm_fraction"].tolist(),
        Cmma_reference_frames=frames, new_DFT_calls=0, physical_parameters_changed=False,
        limitations=["E023 reference choices unchanged: all four tied frames, no G2 coverage-based selection",
            "G2 training observation, not independent prediction or G3 bottleneck selection",
            "common T mass-metric displacement fractions are not affine-origin fractions or energy partitions",
            "atomic displacement and cell strain are separate; an atomic basis does not reconstruct lattice change",
            "source LDA geometric directions are not PBE curvatures, local phonons or new eigenpairs",
            "ordinary residual pass is not a continuous TS or small-barrier-error certificate",
            "raw/transformed author mode arrays are not redistributed"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("case", "source-root", "registration", "variants", "gamma", "output"):
        parser.add_argument("--"+name, type=Path, required=True)
    parser.add_argument("--terminal-audit-sha256", required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("fresh analysis output required")
    result = analyze(args.case, args.source_root, args.registration, args.variants, args.gamma,
                     terminal_audit_sha256=args.terminal_audit_sha256)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, allow_nan=False); stream.write("\n")
    print(json.dumps({k: result[k] for k in ("source_job_id", "replayed_fmax_eV_A", "highest_image_index", "new_DFT_calls")}))


if __name__ == "__main__":
    main()
