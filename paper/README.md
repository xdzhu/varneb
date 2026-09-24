# VARNEB paper workspace

The canonical Computer Physics Communications manuscript is maintained in
[`VARNEB_CPC/`](VARNEB_CPC/).  The former misspelled manuscript directory and
the obsolete root-level `vcneb_CPC` draft have been retired so that the
repository has one manuscript source of truth.

## Layout

- `VARNEB_CPC/varneb_CPC.tex`: current compact CPC manuscript.
- `VARNEB_CPC/varneb.bib`: bibliography used by the manuscript.
- `VARNEB_CPC/figures/`: committed manuscript figures and source-data CSVs.
- `VARNEB_CPC/MANUSCRIPT_EVIDENCE.md`: manuscript-specific evidence checklist.
- `claim_evidence.md`: project-wide claim-to-code-and-result map.
- `references/`: local reference material; intentionally ignored and never
  included in release archives.

## Build

From `paper/VARNEB_CPC/` run:

```text
latexmk -pdf -interaction=nonstopmode -halt-on-error varneb_CPC.tex
```

LaTeX products are local build artifacts.  The manuscript source, bibliography,
editable figures, and source data are versioned.  The article is intentionally
kept near ten typeset pages; new material is added only when it closes a
documented evidence gap.

## Evidence policy

- Every numerical claim must map to a committed source-data file, audit, or
  reproducible calculation record.
- Running or unconverged calculations are described as provisional and do not
  replace accepted production evidence.
- Calculator support and material validation are reported at their demonstrated
  level; a single-image smoke test is not presented as a converged cross-backend
  VCNEB benchmark.
- Large DFT outputs, restart files, private reference packages, and scheduler
  scratch directories remain outside the manuscript package.
