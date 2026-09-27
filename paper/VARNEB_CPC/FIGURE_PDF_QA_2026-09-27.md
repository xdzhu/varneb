# CPC working-PDF and figure QA, 2026-09-27

This is a **working-tree review**, not a submission certificate. The current
TeX project builds successfully with TeX Live/`latexmk` into an isolated
`build_review_20260927/` directory. The rendered manuscript is **12
pages**. Five figure PDFs appear on pages 4, 7, 8, 9, and 10, respectively;
the declarations are on page 11 and the references occupy page 12. The
review build is not substituted for the manuscript package's main PDF.

| Figure | One-sentence conclusion and evidence boundary | Visual/source-data check |
| --- | --- | --- |
| 1, architecture | One manager owns path state and delegates static image evaluations; the listed LAMMPS adapter is not a GaN material-path validation. | Workflow labels and arrows are legible on page 4; schematic only, no quantitative source data. |
| 2, BTO frozen Γ cut | The actual variable-cell T→C path projects onto but leaves a frozen two-soft-mode, fixed-cubic-cell energy cut; this is not a conditional PES or a barrier map. | Rebuilt page 7 with 59 computed static points and seven path values. Panels, regular-weight `(a)` labels, distinct energy zeros, and the off-plane residual are legible. The figure PDF, SVG, PNG, source CSV, and QA record are packaged together. |
| 3, material paths | BTO is monotonic; HfO₂ has an interior maximum whose estimate depends on path protocol. | Four panels and external-literature marker are legible on page 8; path and barrier source CSVs are tracked. Literature values use different calculator contracts and are not error bars. |
| 4, GaN backends | Five complete 45.7-GPa paths share the dominant barrier topology; forward and reverse bars use different endpoint baselines. | Three panels are legible on page 9; all five backends have 29 rows in the tracked source CSV. The Qian curve is an approximate digitization, not a sixth computed backend. |
| 5, local GaN coupling | A negative direction appears in the combined atom–strain Hessian even though its frozen blocks are positive. | Two panels and the signed-curvature zero line are legible on page 10; tracked CSV and audit accompany the figure. The 1000-eV joint Hessian and 600-eV endpoint Γ projections are compared geometrically, not as one energy protocol or decomposed barrier. |

The plotted SVGs retain text elements; the GaN figure PDFs have embedded
TrueType text rather than text converted entirely to outlines. Top/right axes,
regular-weight panel labels, and framed legends are present. No obvious plot
line/legend collisions or clipped labels were observed in the reviewed PNGs
or rendered manuscript pages. This visual check does **not** establish the
physical accuracy of an individual DFT point; those claims require the
separate raw-output audits.
The BTO figure uses display-only interpolation. Final-grid interior
leave-one-out error reaches 5.14 meV/BTO; prior independent coarser-grid
center and edge holdouts reached 32.57 and 14.83 meV/BTO. None of those
checks certifies unsampled extrema. The previously generated pure mode-norm
figure remains archived but is not a main-text figure in this build.

The first review render exposed one print-size defect: Table 3's long
`BlockFIRE (0.02/image)` labels touched the step counts. The table now uses
short labels and defines their per-image step cap in its caption. A second
successful build and page-7 render verified that the columns no longer touch.
Two math-containing subsection bookmarks were given plain-text PDF strings;
the latest build has no `Token not allowed` bookmark warnings. One 1.9-pt
overfull output box remains and requires a final print-size review.
The working tree's architecture PDF/SVG were already modified before this
review; they were neither staged nor treated as part of the new BTO figure.

Remaining submission gates:

- Trim or justify the 12-page draft against the approximately 10-page target
  after the BTO holdout and final figure set are frozen. The reference-only
  last page makes layout/float economy worth reviewing, but scientific methods
  and limitations must not be removed merely to save space.
- Confirm author order, affiliations, funding, CRediT roles, declarations,
  and the CPC Library entry; page 1 and page 11 still contain placeholders.
- Freeze hashes for every final figure PDF/CSV and reconcile all captions,
  backend parameters, normalization, and bibliography against committed
  evidence. Several concurrent figure/manuscript changes remain uncommitted
  in the shared worktree and are not silently included in a release claim.
