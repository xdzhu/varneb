"""No-DFT regression test for the Phonopy Gamma-eigenvector export CLI."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_exporter_declares_a_no_dft_phonopy_qpoint_contract() -> None:
    text = (ROOT / "scripts" / "export_phonopy_gamma_eigenpairs.py").read_text(encoding="utf-8")
    for fragment in (
        "from phonopy import load",
        "phonopy_gamma_eigenpairs(phonon)",
        "force_sets_filename=str(force_sets_path)",
        "is_compact_fc=False",
        "with_eigenvectors=True",
        "phonopy_yaml_sha256",
        "force_sets_sha256",
    ):
        assert fragment in text
