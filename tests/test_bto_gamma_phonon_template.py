"""Regression guards for the BTO Gamma phonon seed workflow."""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_bto_gamma_phonon_template_is_static_and_opt_in() -> None:
    text = (ROOT / "cluster" / "hf_batio3_gamma_phonon_abacus.slurm").read_text(encoding="utf-8")
    for fragment in (
        "#SBATCH --ntasks=32",
        "--abacus -d --pm --dim \"1 1 1\"",
        "calculation scf",
        "cal_force 1",
        "cal_stress 1",
        "symmetry 0",
        "if [[ \"${RUN_DFT:-0}\" != 1 ]]",
        "RUN_DFT=0 so no ABACUS SCF was launched",
        "phonon_preflight.json",
        "'endpoint': __import__('os').environ['ENDPOINT']",
        "--abacus -f disp-[0-9][0-9][0-9]/OUT.ABACUS/running_scf.log",
        "force_sets_filename='FORCE_SETS',",
        "calculator='abacus',",
        "is_compact_fc=False",
        "write_FORCE_CONSTANTS(phonon.force_constants)",
        "transpose(0, 2, 1, 3)",
        "gamma_force_constants.npz",
        "phonopy_gamma_eigenpairs.npz",
        "phonopy_gamma_eigenpairs(phonon)",
        "save_phonopy_gamma_eigenpairs",
        "force_constant_provenance.json",
    ):
        assert fragment in text
