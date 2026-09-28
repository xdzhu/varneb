"""Regression checks for the read-only CP2K raw/cache audit."""

from pathlib import Path
import os

from ase.io import read
import numpy as np

from scripts.audit_cp2k_image_output_cache import compare
from scripts.audit_cp2k_final_chain_raw_energy import (
    ENERGY_RE,
    cache_digest,
    last_stress_gpa,
)


ROOT = Path(__file__).resolve().parents[1]


def test_exact_final_peak_cache_key_matches_archived_chain():
    trajectory = (ROOT / "paper/VARNEB_CPC/evidence"
                  / "gan_45p7_final_chains_20260927"
                  / "gan_cp2k_45p7_final_chain.traj")
    peak = read(trajectory, index=15)
    assert cache_digest(15, peak) == (
        "cc482fbb3276f4d82daccec61e40f26c4c6d87a17e4566c10f7016afadef530d"
    )


def test_last_stress_block_and_force_eval_energy_parsing():
    text = """
 ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]: -169.123456789000
 STRESS|                        x                   y                   z
 STRESS|      x        4.100000E+01   0.000000E+00   0.000000E+00
 STRESS|      y        0.000000E+00   4.200000E+01   0.000000E+00
 STRESS|      z        0.000000E+00   0.000000E+00   4.300000E+01
 ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]: -169.056699139691062
 STRESS|                        x                   y                   z
 STRESS|      x        4.68941641381E+01  -3.79007101902E-05  -2.37509620080E-05
 STRESS|      y       -3.79007101902E-05   4.47196007532E+01  -2.61781073970E-03
 STRESS|      z       -2.37509620080E-05  -2.61781073970E-03   4.53809748806E+01
"""
    assert len(ENERGY_RE.findall(text)) == 2
    stress = last_stress_gpa(text)
    assert stress.shape == (3, 3)
    assert abs(stress[0, 0] - 46.8941641381) < 1e-11
    assert abs(stress[1, 2] + 0.00261781073970) < 1e-13


def test_cumulative_output_cache_audit_detects_unit_factor(tmp_path):
    output = tmp_path / "cp2k.out"
    output.write_text(
        " ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]: -169.100000000000000\n"
        " ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]: -169.050000000000000\n",
        encoding="utf-8",
    )
    cache_dir = tmp_path / "image_cache"
    cache_dir.mkdir()
    for index, raw_energy in enumerate((-169.1, -169.05)):
        path = cache_dir / f"image_0015_{index:064x}.npz"
        np.savez_compressed(path, energy=raw_energy * 27.21138385655639)
        timestamp_ns = (index + 1) * 1_000_000_000
        os.utime(path, ns=(timestamp_ns, timestamp_ns))
    report = compare(output, cache_dir, 15, 1e-8)
    assert report["n_output_records"] == 2
    assert report["n_unique_cache_entries"] == 2
    assert report["n_matching_cache_entries_with_ase_constant"] == 0
    assert abs(report["effective_order_paired_hartree_to_ev"]
               - 27.21138385655639) < 1e-12
    assert report["max_absolute_energy_residual_after_factor_ev"] < 1e-10
