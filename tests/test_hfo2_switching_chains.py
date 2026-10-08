from pathlib import Path
import json

from ase.io import read
import numpy as np
import pytest

from scripts.prepare_hfo2_switching_chains import ordered_seed, prepare
from vcneb import endpoint_structure_record, validate_path_geometry


def test_real_ordered_variants_are_distinct_safe_seeds_and_exact_endpoints():
    root = Path(__file__).resolve().parents[1] / "benchmarks/hfo2_channels/20261008/reference_variants"
    initial = read(root / "PO.vasp")
    chains = []
    for name in ("PO_minus_T_preserving_inversion", "PO_minus_T_reversing_inversion"):
        final = read(root / f"{name}.vasp")
        chain = ordered_seed(initial, final)
        assert len(chain) == 9
        # VCNEB optimizes unwrapped coordinates. Replacing the MIC-lifted
        # final image by its wrapped input produces a false final long segment.
        segments = validate_path_geometry(chain)["segment_lengths_A"]
        assert np.allclose(segments, segments[0], atol=1e-10, rtol=1e-10)
        fractional_steps = np.diff(np.array([
            image.get_scaled_positions(wrap=False) for image in chain
        ]), axis=0)
        assert np.allclose(fractional_steps, fractional_steps[0], atol=1e-12, rtol=0)
        final_shift = chain[-1].get_scaled_positions(wrap=False) - final.get_scaled_positions(wrap=False)
        assert np.allclose(final_shift, np.rint(final_shift), atol=1e-12, rtol=0)
        for expected, actual in ((initial, chain[0]), (final, chain[-1])):
            assert endpoint_structure_record(expected)["sha256"] == endpoint_structure_record(actual)["sha256"]
        for image in chain:
            assert image.get_chemical_symbols() == initial.get_chemical_symbols()
            assert image.calc is None
            distances = image.get_all_distances(mic=True)
            np.fill_diagonal(distances, np.inf)
            assert distances.min() > 2.0
        chains.append(chain)
    assert np.linalg.norm(chains[0][4].positions - chains[1][4].positions) > .5
    changed = initial.copy()
    changed.cell[0, 0] += .01
    with pytest.raises(ValueError, match="same-cell"):
        ordered_seed(initial, changed)


def test_no_switching_seed_before_real_electronic_gate(tmp_path):
    polar = tmp_path / "polar"
    polar.mkdir()
    (polar / "manifest.json").write_text("{}")
    (polar / "summary.json").write_text(json.dumps({"status": "review_required"}))
    with pytest.raises(ValueError, match="not passed"):
        prepare(tmp_path, polar, tmp_path / "not_created")
    assert not (tmp_path / "not_created").exists()


def test_prepared_hf_seed_evidence_keeps_byte_hashes_and_continuous_lifts():
    from scripts.audit_hfo2_static_replica import sha256
    from vcneb import validate_periodic_path_lift
    root = Path(__file__).resolve().parents[1] / "benchmarks/hfo2_channels/20261008"
    preflight = json.loads((root / "switching_preflight_r11.json").read_text())
    assert preflight["status"] == "passed" and not preflight["DFT_executed"]
    for record in preflight["reports"]:
        folder = root / "switching_seeds" / record["variant"]
        manifest = json.loads((folder / "manifest.json").read_text())
        assert sha256(folder / "seed.traj") == record["seed_traj_sha256"]
        assert sha256(root / "polarization_summary.json") == manifest["polarization_summary_sha256"]
        assert sha256(root / "polarization_manifest_r10.json") == manifest["polarization_manifest_sha256"]
        images = read(folder / "seed.traj", index=":")
        validate_periodic_path_lift(images)
        assert len(images) == 9 and record["n_active_images"] == 7
        assert all(a["new_SCFs"] == 0 for a in record["endpoint_cache_audits"])
