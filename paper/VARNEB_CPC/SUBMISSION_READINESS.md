# VARNEB CPC submission gate

Current scope: one manuscript with no page limit and no separate supplement.
The 2026-09-30 build has eight main-text figures and one same-document
appendix figure. This gate supersedes the layout targets in the dated
`STAGE_GATE_2026-09-27.md`; it does not supersede raw-data audits.

## CPC article route

The present Program Summary and installable-software emphasis make this a
candidate *Computer Programs in Physics* (CPiP) submission. The
[publisher's journal description](https://shop.elsevier.com/journals/computer-physics-communications/0010-4655)
distinguishes CPiP from *Computational Physics Papers*: CPiP software is
archived in the CPC Program Library on Mendeley Data and must have an approved
open-source licence. The repository declares GPL-3.0-or-later in
`pyproject.toml` and includes the GPLv3 licence text, but the CPC archive has
not yet been deposited and no library link exists. Do not replace that field
with a guessed DOI. The same journal description asks authors to articulate
novelty and physics significance for a general audience; the manuscript must
therefore lead with the portable execution/evidence/mode-analysis capability
and its demonstrated GaN/BTO use, not suggest a new VCNEB physical formalism.
The authors should confirm the CPiP route and archive requirements during
submission rather than silently treating a GitHub or PyPI link as the CPC
Library accession.

## Claims the present evidence can support

| Manuscript claim | Evidence to check | Allowed conclusion and boundary |
| --- | --- | --- |
| Established NEB extended to variable cell through a backend-independent controller | `varneb_CPC.tex` Secs. 2–3; `tests/check_vcneb_forces.py`; fixed-cell and finite-difference regressions | A tested implementation and execution contract, **not** a new VCNEB physical formalism. |
| BTO/ABACUS barrierless T→C path and cubic Γ decomposition | Figs. 2–3 source CSV/QA, `outputs/batio3_t_to_c_pbe100_dzp10au/` provenance, `bto_qy0_27_node_sheet_20260929.json`, `bto_transverse_soft_frozen289_20260929.json` | Under PBE, 100 Ry, 10-au DZP, and a 4×4×4 electronic mesh, the seven-total-image path is monotonic; the 1×1×1 Γ basis explains atomic geometry. The 27-point `Q_y=0` restricted variable-cell sheet and 289-point frozen-cubic two-soft-mode cut are different objects. Their prospectively sampled increments pass 1.131- and 0.488-meV/BTO maximum-error checks, respectively. Neither is an unrestricted conditional PES, finite-temperature free energy, or activation barrier. |
| GaN 45.7-GPa B4→B1 multi-backend barrier topology | Fig. 4 source data, `evidence/gan_45p7_multibackend_vcneb_20260924.json`, `evidence/gan_45p7_full_image_input_contract_audit_20260929.md`, `benchmarks/numerical_integrity/gan_cp2k_final_chain_raw_energy_stress_20260928.json` | Five 29-total-image paths satisfy the common 0.10-eV/Å generalized-force criterion and share a dominant peak. The final visible inputs have consistent per-backend settings, but this is not a historical-call or raw-output audit. Cross-backend energies use separate endpoint baselines and calculator contracts. The 45.7-GPa pressure is adopted from Qian, not re-established as each PBE backend's coexistence pressure; Qian's PW91/ultrasoft/force-norm protocol differs. CP2K original-text atomic forces and final run-end markers are unavailable; cache-chain forces and raw energy/stress reconciliation do not replace them. |
| GaN atom–strain and mode interpretation | Figs. 5–8 source CSV/QA, `gan_600eV_joint_gamma_bridge_20260928.json`, `gan_600eV_ts_2d_17x17_refinement_20260930.json`, `gan_600eV_atomic_tube_q9_20260929.json` | At the original VASP/PBE 600-eV, 45.7-GPa contract, the local joint Hessian has a negative direction. The independently sampled 289-point local atom–strain cut and 162-point central path-transverse cut support a bounded mechanism interpretation. The prior 9×9 surface predicts 208 new joint-cut DFT points within 0.00698 meV/GaN. Two measured transverse points are about 0.2 meV/GaN below their frozen-cell centerline references; they are not a fully relaxed lower MEP. The energy–force/stress mismatch and incomplete stationary-point test preclude a strict variable-cell TS certificate. The central cut is not a whole-path or relaxed two-mode surface. |
| Case-specific convergence acceleration | Table 2, `ACCELERATION_EVIDENCE_AUDIT_2026-09-27.md`, Slurm launch audit | Matched-start ABACUS chains use the same first-crossing threshold of 0.10 eV/Å; BTO/HfO₂ counted-launch reductions are 78.5%/38.2%. Negative transfers and path differences remain visible. Do not advertise universal speedup or identical saddle basins. The pre-DFT feasibility gate is a separate safety mechanism; this ablation did not isolate its contribution to launch savings. |
| Installable and usable source package | `docs/RELEASE_NOTES_v0.0.2.md` hosted-artifact smoke; `docs/USER_MANUAL.md`; `tests/test_cli_run.py`, `tests/test_config.py` | PyPI 0.0.2 passes an isolated CLI smoke but predates the later QE/CP2K/ABINIT evidence. The current source is now numbered 0.0.3, which has separately passed isolated local wheel/sdist build, Twine checks, fresh wheel install, and the sdist analytic toy. No 0.0.3 tag or publication exists; neither artifact smoke certifies a user's DFT executable. |

`MANUSCRIPT_EVIDENCE.md` and `FIGURE_LOGIC_AND_STYLE.md` are the current
figure-to-source maps. A pretty interpolated pixel is never an additional
DFT observation. The claim limits above must remain in the abstract, results,
captions, and conclusions when wording changes. The direct manuscript-number
regression `tests/test_cpc_numeric_claims.py` recomputes headline GaN barriers,
BTO blind errors and monotonic endpoint rise, Slurm launch reductions, and
the HfO$_2$ per-cell conversion from committed source tables.

## Clean-source QA snapshot, 2026-09-30

The pushed `037ace6` source exposed a packaging gap that the shared working
tree hid: the four-landscape rebuild used three small BTO numeric provenance
JSON files under ignored `outputs/`. Its first independent Git-archive test
had one failure (602 passed, one skipped); no DFT or figure-audit result was
invalidated. Commit `8ba4971` added exactly those three source inputs without
changing their bytes or any ABACUS calculation parameter. From a fresh
`git archive` of `8ba4971`, the suite passes **603 tests, one skipped**;
the one-command rebuild regenerates all four main BTO/GaN landscapes with
source CSVs and PNGs byte-identical to the archived artwork, including the
new GaN 17×17 cut. The archive's LaTeX build passes and yields 15 A4 pages.
The shared working-tree suite after the GaN-figure update passes **622 tests,
one skipped**; its extra tests include unrelated untracked local work and
cannot substitute for the clean-source check. These checks validate the
packaged data-to-figure pipeline, not a new DFT or index-one TS certificate.

The bounded VASP-transfer wording and evidence-map update passed the shared
working-tree suite (**624 passed, one skipped**). It initially built as 16 A4
pages with a sparse two-column conclusions page and a separate nearly empty
appendix-introduction page. A forced page break did not help and was reverted.
The subsequent local layout edit balances the last two-column page and places
the appendix text beside its figure in one-column backmatter; TeX Live now
builds **15 A4 pages**, and pages 13--15 were visually inspected. Two existing
minor overfull boxes remain (1.9 and 0.58 pt); there are no unresolved
references or citations in the build log. Author metadata and the final-source
clean build remain open gates; this working-tree build includes an unrelated
locally modified architecture figure and must not be called a frozen
submission PDF.

An independent `git archive` of the pushed `1176eb3` revision now closes the
source-only check for this draft stage: **605 tests passed, one skipped** from
the extracted tree; TeX Live independently compiled the archived source to
**16 A4 pages**. Its one-command four-landscape rebuild reproduced all four
archived source CSVs and PNGs byte for byte, with the expected 27, 289, 289,
and 162 measured-node counts. The locally modified architecture PDF differs
in file bytes from `1176eb3`, but both versions rasterize to an identical
144-dpi PNG at the same page size; that visual equivalence does not stage or
approve the local edit. This is a reproducibility check for a provisional
revision, not the final CPC Library deposit or a replacement for author review.

The pushed `c754714` layout revision independently passes **605 tests, one
skipped** and builds a **15-page** PDF from a fresh Git archive; its rendered
pages 13--14 match the working-tree inspection pixels. A temporary virtual
environment created inside that archive ran the README's documented editable
install with normal build isolation (`pip install -e ".[plot]"`), reporting
VARNEB 0.0.3. The analytic toy returned a 0.250004-eV barrier. The bundled
ASE/EMT JSON then passed `validate-config`, `prepare`, and `run --execute`:
three total images, zero barrier, and final force
`1.52e-18 eV/Å`; a second execution was correctly refused. `backends --json`,
`optimizers --json`, `doctor`, and the `v0.0.3` metadata check also succeeded.
The virtual environment inherited ASE/NumPy from the host, so this is a
clean-source CLI smoke, not a fully isolated dependency-resolution test or a
DFT-backend certification. A deliberately non-documentary
`--no-build-isolation` attempt failed against the host's old setuptools;
the documented build-isolated install succeeded without changing package
metadata. No DFT executable was present on this local `PATH`.

On 2026-09-30, the two GaN two-dimensional figures received explicit,
non-obscuring star-marker legends. Fig. 5's blue star is a one-step-refined
image-15 static and Fig. 8's orange star is the original path's highest
discrete image 15; neither is an endpoint or a strict TS certificate. Their
source CSV bytes and DFT counts are unchanged. The shared-tree suite now
passes **626 tests, one skipped**; an independent four-landscape rebuild
reproduces all archived source CSVs and PNGs byte for byte with the expected
27/289/289/162 measured nodes. The revised TeX compiles to **15 A4 pages**,
and the affected figure pages 10 and 12 were rendered and visually checked.
These are figure/claim-clarity and reproducibility checks, not new DFT.

A subsequent hash-bound read of the original VASP production image-0/28
static OUTCARs, now archived without POTCAR, corrects an endpoint-provenance
ambiguity: original B4/B1 raw-stress residuals are **2.912/0.239 kbar** at
45.7 GPa, whereas 1.743/2.528 kbar belong to later, structurally close
signed basin-return relaxations. Both original endpoint forces are below
0.02 eV/Å and their raw energies match the final chain exactly. The B4
production endpoint misses the separately adopted 2-kbar diagnostic gate;
the present paper reports this rather than claiming endpoint stress
certification or silently replacing its endpoint.

On 2026-10-01, a bounded same-electronic-contract B4 endpoint check resolved
this narrow miss: the first native relaxation (27812043) stopped immediately
at `EDIFFG=-0.02 eV/A` and retained 2.912 kbar; its stricter isolated
continuation (27812112, `EDIFFG=-0.005 eV/A`) took five evaluated frames and
reached 0.059 kbar with 0.000151 eV/Å maximum atomic force. Relative atomic
displacement from the original B4 is at most 0.00225 Å and the native
enthalpy change is -0.036 meV/GaN. The jobs did not modify the production
chain, its endpoint cache, the figures, or the reported barrier. They show
small numerical endpoint slack, not a different basin or a strict TS
certificate. The local-cut star is an image-15 reference, not an endpoint.

A further 12-case same-600-eV normal-strain step-size audit (Slurm array
`27812251`) finds that the image-15 energy–stress residual persists at
0.0205–0.0238 eV/Å across 0.02, 0.01, and 0.005 Å steps. This rules out a
simple large-step truncation explanation over that tested range, not every
possible numerical cause. Individual k-point plane-wave counts change in
355–361/384 k points even for the 0.005-Å pairs, consistent with but not
uniquely proving a finite-basis contribution. The manuscript now reports
this evidence and keeps the strict TS claim explicitly unproven; the
enthalpy landscape and original barrier are unchanged.

The GaN Fig. 8 negative off-center points now also have a near-center
energy–force consistency check: image-19/20 ±0.0125-Å secants differ from
projected path-force slopes by only 0.0002893/0.0000017 eV/Å. The two
centers are 117.54/158.50 meV/GaN below the discrete peak, and their
0.237/0.205-meV/GaN frozen lowering does not change that peak. This
strengthens the sampled-landscape interpretation while leaving an alternative
fully relaxed MEP unproven; it is not a reason to rerun the full chain.

## Historical technical QA snapshot, 2026-09-29

On the shared working tree after the 289/27/162-point figure update,
`python -m pytest -q --disable-warnings` completed with **618 passed, one
skipped**, and 209 suppressed warnings (mostly ASE/spglib/Phonopy
deprecations). The first unscoped attempt collected duplicate test modules
from retained `tmp/` source snapshots and failed during collection; the
repository now declares `testpaths = ["tests"]`, and both the explicitly
scoped and plain-root entry points pass. Figure regression tests have been
updated to check the current main-text CSVs against the 27/289/162-node
audit records rather than merely expecting old filenames.
An independently extracted clean Git archive containing the 27/289/162-point
figures passed **599 tests, one skipped**, and compiled the 15-page manuscript
without missing figures or citations. This checks committed source rather
than relying on the shared working tree. Its smaller test count reflects
untracked tests in that tree, not failures. The first archive audit exposed
that Git normalized three CRLF raw-audit JSON files to LF, breaking their
figure-QA SHA-256 links; exact-file `.gitattributes` rules now preserve those
bytes. Current regressions check the report, plotter, and source-data hashes.
Increasing the
permitted double-column top-float occupancy then packed GaN Figs. 5--8
beside their discussion across pages 9--11; a fresh working-tree TeX Live
build after the new figures is **15 A4 pages**. Pages 6, 7, and 10–13 were
rendered and checked for clipping;
the remaining appendix/reference-page whitespace and one 1.9-pt overfull
box are layout polish, not missing calculations. Moving the main-text float
barrier past the conclusion instead placed figures after their discussion,
so that earlier trial was reverted. This is not a final-release test: the
shared working tree contains unrelated modified figures and untracked tests,
and author/submission metadata will still change the document.
The documented ASE/EMT quickstart also passed in a copied isolated directory:
the public CLI now emits one concise line while preserving the complete
`vcneb_summary.json`; `--full-summary` retains the old verbose option.
The exact `519a697` Git-source archive independently ran the README analytic
toy (0.250004-eV barrier) and the bundled ASE/EMT JSON quickstart through
`python -m vcneb validate-config`, `prepare`, and `run --execute`. The latter
converged in its one-step interface check with a zero barrier and
`fmax=1.52e-18 eV/Å`; it did not launch DFT. The installed
`varneb` console entry point is checked separately by the wheel smoke.
After distinguishing the post-0.0.2 source as version 0.0.3, a new local
isolated build produced wheel and sdist artifacts that passed `twine check`.
The wheel installed in a fresh Python 3.10 virtual environment, reported
version 0.0.3, and exposed the expected CLI; the unpacked sdist toy returned
0.250004 eV. `tests/check_release_metadata.py --tag v0.0.3` and the full
working-tree suite (616 passed, one skipped) passed. This is an unpublished
candidate, not a CPC Library deposit or PyPI release.

A targeted manuscript artifact check found all nine referenced PDF figures
and nine corresponding scientific source CSVs present and Git-tracked. The
TeX source has 27 labels and 15 unique references with no missing target;
19 unique citation keys resolve against 20 BibTeX entries. This is a
presence/linkage check, not a raw-DFT audit. Figure 1's architecture PDF and
SVG are locally modified in the shared working tree (the SVG diff is limited
to generation metadata and clip IDs); they were not staged or pushed here.
The final frozen manuscript build must use a reviewed, committed version of
that figure so its bytes match the archived source revision.

The locally built 0.0.3 PyPI source distribution is an **installable package**,
not the complete CPC research archive. Its 51 tar entries include 31 `vcneb`
Python modules, the user manual, and one analytic toy example, but no tests,
material-case inputs, or manuscript source data; this follows `MANIFEST.in`.
The CPiP program deposit should therefore be prepared from a pinned final
repository revision (or an equivalently complete, reviewed source snapshot),
with the package, tests, representative runnable examples, and the figure/data
provenance index together. Keep the copyrighted DFT executables, licensed
potentials, and nonredistributable raw inputs out of that public snapshot;
verify its contents and licence before deposit. The PyPI sdist smoke proves
installation, not reproduction of the material-level paper figures.

An independent `git archive` extraction exposed a repository-only failure
that the shared working-tree test had hidden: the first clean snapshot had
12 failures (585 passed, one skipped). Raw-byte SHA-256 links broke when
archive line-ending conversion changed hash-bound CSV, JSON, Python, and
Phonopy text files; the BTO restricted-sheet test also needed seven small
provenance/structure/trajectory inputs present locally but omitted from Git by
the local `outputs/` exclusion. A staged-source candidate now preserves those
raw bytes through explicit `.gitattributes` rules, tracks only the seven
required BTO files, and passes its independent archive suite (597 passed,
one skipped). That count was lower than the corresponding shared-tree
snapshot because it included untracked tests from another work stream. The
newer clean-archive check above supersedes this historical count. Neither
archived source test is a fresh DFT audit or a final CPC/Mendeley deposit.

## Remaining gates for this CPC submission

1. Replace author-order, affiliation, corresponding-author, CRediT,
   funding/computing-acknowledgement, and competing-interest placeholders
   with author-approved statements. Confirm the CPiP article route and
   resolve the Program Summary's pending CPC Library field through the actual
   CPC/Mendeley submission workflow.
2. The four quantitative external comparators have a primary-page read-through
   in `FINAL_PRIMARY_SOURCE_READTHROUGH_2026-09-29.md`; the BTO PBEsol
   distortion versus NEB distinction and Qian/Sheppard/Liu protocol limits
   are retained. Before submission, obtain an author-level final read of
   captions, pressure and formula-unit normalization, the CP2K and TS limits,
   general-method citations, and the Data availability statement. The DOI
   metadata and corrected author records are separately documented in
   `CITATION_AUDIT_2026-09-28.md`.
3. On the final source revision, rerun the full test suite, manuscript claim
   check, LaTeX build, and rendered-page inspection. Review the late-figure
   spacing without moving figures past the conclusions. Record failures
   rather than treating an old green result as proof for changed files.
4. Pin the exact final code revision in the CPC program deposit and data
   availability record. PyPI 0.0.2 is an earlier snapshot; 0.0.3 is locally
   built but unpublished. Publish a reviewed 0.0.3 release only when the
   source and submission package are frozen, or cite the final repository
   commit explicitly. Prepare and inspect a complete CPC source snapshot
   separately from the deliberately minimal PyPI sdist. Do not identify
   0.0.2 as the code state used for the later GaN evidence.

No new material calculation is a prerequisite for the **bounded claims
actually made in this draft**, unless a final source audit finds an error.
The unrestricted BTO lower envelope is obstructed at cubic C by the free
third soft Γ direction; strict GaN index-one TS certification, HfO₂ multi-mode
maps, bilayer hBN, and mode-guided release-and-refine remain follow-on
research, not retroactive claims of this CPC paper. Do not adjust BTO's
100-Ry/10-au-DZP contract or GaN's 600-eV contract to improve a plot.
