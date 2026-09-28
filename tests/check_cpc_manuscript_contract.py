"""Guard the material-claim boundaries in the CPC manuscript source."""

from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "paper" / "VARNEB_CPC" / "varneb_CPC.tex"
SUPPLEMENT = ROOT / "paper" / "VARNEB_CPC" / "varneb_CPC_supplement.tex"
BIBLIOGRAPHY = ROOT / "paper" / "VARNEB_CPC" / "varneb.bib"


def require(text: str, fragment: str, *, label: str) -> None:
    # LaTeX prose is manually wrapped; a line break is not a scientific change.
    normalized_text = re.sub(r"\s+", " ", text)
    normalized_fragment = re.sub(r"\s+", " ", fragment)
    if normalized_fragment not in normalized_text:
        raise SystemExit(f"missing {label}: {fragment!r}")


def main() -> None:
    manuscript = MANUSCRIPT.read_text(encoding="utf-8")
    supplement = SUPPLEMENT.read_text(encoding="utf-8")
    bibliography = BIBLIOGRAPHY.read_text(encoding="utf-8")

    # These are scientific reporting boundaries, not stylistic preferences.
    for fragment, label in (
        ("without changing the underlying VCNEB formalism", "algorithm-novelty boundary"),
        ("Material-level agreement requires a recorded", "backend evidence boundary"),
        ("same dominant barrier topology", "cross-backend claim boundary"),
        ("there is no interior barrier, so CI is withheld", "BTO topology and CI policy"),
        ("preclude a matched algorithm benchmark", "BTO literature qualifier"),
        ("0.1291722", "HfO2 CI barrier evidence"),
        ("reported 32 meV per formula unit", "HfO2 literature comparison"),
        ("78.5\\% and 38.2\\%", "bounded acceleration result"),
        ("(b) Forward $H_{\\rm peak}-H_{\\rm B4}$", "forward GaN barrier panel"),
        ("(c) Reverse $H_{\\rm peak}-H_{\\rm B1}$", "reverse GaN barrier panel"),
        ("At the 2-kbar stress gate", "GaN basin stress gate"),
        ("B1 narrowly misses (2.528 kbar)", "GaN B1 stress limitation"),
        ("not a stress-certified TS", "GaN TS claim limit"),
        ("81 measured enthalpies", "GaN local DFT point count"),
        ("or whole-path surface or a certified stationary transition state",
         "GaN local TS limit"),
        ("Supplementary Fig.~S4", "BTO frozen-grid supplement boundary"),
    ):
        require(manuscript, fragment, label=label)

    require(supplement, "vcneb_material_validation.pdf", label="material validation supplement")
    require(supplement, "bto_frozen_soft_mode_81_20260928.pdf",
            label="audited BTO grid supplement")
    require(manuscript, "\\cite{BTOReaxFF2019}", label="BTO citation key")
    require(bibliography, "@article{BTOReaxFF2019,", label="BTO bibliography record")
    require(bibliography, "10.1039/C9CP02955A", label="BTO DOI")
    print("cpc_manuscript_contract=ok")


if __name__ == "__main__":
    main()
