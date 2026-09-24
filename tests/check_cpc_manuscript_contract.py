"""Guard the material-claim boundaries in the CPC manuscript source."""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "paper" / "VARNEB_CPC" / "varneb_CPC.tex"
BIBLIOGRAPHY = ROOT / "paper" / "VARNEB_CPC" / "varneb.bib"


def require(text: str, fragment: str, *, label: str) -> None:
    if fragment not in text:
        raise SystemExit(f"missing {label}: {fragment!r}")


def main() -> None:
    manuscript = MANUSCRIPT.read_text(encoding="utf-8")
    bibliography = BIBLIOGRAPHY.read_text(encoding="utf-8")

    # These are scientific reporting boundaries, not stylistic preferences.
    for fragment, label in (
        ("without changing the underlying VCNEB formalism", "algorithm-novelty boundary"),
        ("Material-level agreement requires a recorded", "backend evidence boundary"),
        ("same dominant\nbarrier topology", "cross-backend claim boundary"),
        ("no interior\nbarrier", "BTO topology conclusion"),
        ("CI is correctly withheld", "BTO CI policy"),
        ("not a strict algorithm benchmark", "BTO literature qualifier"),
        ("0.1291722", "HfO2 CI barrier evidence"),
        ("reported 32 meV per formula unit", "HfO2 literature comparison"),
        ("figures/vcneb_material_validation.pdf", "reproducible figure asset"),
        ("78.5\\% and 38.2\\%", "bounded acceleration result"),
    ):
        require(manuscript, fragment, label=label)

    require(manuscript, "\\cite{BTOReaxFF2019}", label="BTO citation key")
    require(bibliography, "@article{BTOReaxFF2019,", label="BTO bibliography record")
    require(bibliography, "10.1039/C9CP02955A", label="BTO DOI")
    print("cpc_manuscript_contract=ok")


if __name__ == "__main__":
    main()
