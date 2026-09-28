"""Archived BTO soft triplet has a free unstable axis outside the 2D plane."""

from pathlib import Path

import numpy as np
import pytest

from scripts.audit_bto_soft_triplet_conditional_endpoint import audit, anchor_projection


ROOT = Path(__file__).resolve().parents[1]
BTO = ROOT / "outputs" / "batio3_t_to_c_pbe100_dzp10au"
FIGURES = ROOT / "paper" / "VARNEB_CPC" / "figures"


def test_archived_cubic_endpoint_obstruction_and_measured_branches():
    report = audit(
        BTO / "bto_cubic_phonopy_gamma_eigenpairs.npz",
        BTO / "bto_cubic_gamma_force_constants.npz",
        BTO / "bto_cubic_phonopy_gamma_eigenpairs_provenance.json",
        BTO / "bto_cubic_gamma_FORCE_SETS",
        FIGURES / "bto_conditional_four_point_stage_2026-09-27_source_data.csv",
        FIGURES / "bto_frozen_soft_mode_landscape_source_data.csv",
    )
    assert len(report["imaginary_triplet_frequencies_THz_signed"]) == 3
    assert max(report["imaginary_triplet_frequencies_THz_signed"]) < 0
    assert report["anchor_projection"]["omitted_y_fraction_outside_fixed_zx"] == pytest.approx(1)
    assert report["archived_T_projection_Qz_Qx_sqrt_amu_A"][0] == pytest.approx(1.20432987766)
    assert report["archived_C_projection_Qz_Qx_sqrt_amu_A"] == pytest.approx([0, 0], abs=1e-10)
    assert report["branch_lowering_range_meV_per_BTO"] == pytest.approx(
        [28.30807676900804, 57.46661343391679]
    )


def test_omitted_anchor_must_be_independent_of_fixed_plane():
    # An artificial triplet with no Ba_y/Ti_y content must fail the audit.
    permutation = [0, 2, 6] + [i for i in range(15) if i not in (0, 2, 6)]
    vectors = np.eye(15)[:, permutation]
    with pytest.raises(ValueError, match="not independent"):
        anchor_projection(vectors, np.ones(5))
