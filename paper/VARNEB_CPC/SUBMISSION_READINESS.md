# VARNEB CPC submission gate

Current scope: one manuscript with no page limit and no separate supplement.
The 2026-09-28 build has eight main-text figures and one same-document
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
| BTO/ABACUS barrierless T→C path and cubic Γ decomposition | Fig. 2/3 source data, `outputs/batio3_t_to_c_pbe100_dzp10au/` provenance, `bto_qy0_even_mode_two_holdout_gate_20260928.json` | Under PBE, 100 Ry, 10-au DZP, and a 4×4×4 electronic mesh, the seven-total-image path is monotonic; the 1×1×1 Γ basis explains atomic geometry. The 81-point frozen cut and nine-training/two-prospective-point `Q_y=0` restricted sheet are different objects. Neither is an unrestricted conditional PES, finite-temperature free energy, or activation barrier. |
| GaN 45.7-GPa B4→B1 multi-backend barrier topology | Fig. 4 source data, `evidence/gan_45p7_multibackend_vcneb_20260924.json`, `evidence/gan_45p7_full_image_input_contract_audit_20260929.md`, `benchmarks/numerical_integrity/gan_cp2k_final_chain_raw_energy_stress_20260928.json` | Five 29-total-image paths satisfy the common 0.10-eV/Å generalized-force criterion and share a dominant peak. The final visible inputs have consistent per-backend settings, but this is not a historical-call or raw-output audit. Cross-backend energies use separate endpoint baselines and calculator contracts. The 45.7-GPa pressure is adopted from Qian, not re-established as each PBE backend's coexistence pressure; Qian's PW91/ultrasoft/force-norm protocol differs. CP2K original-text atomic forces and final run-end markers are unavailable; cache-chain forces and raw energy/stress reconciliation do not replace them. |
| GaN atom–strain and mode interpretation | Figs. 5–8 source CSV/QA, `gan_600eV_joint_gamma_bridge_20260928.json`, `gan_600eV_ts_2d_9x9_refinement_20260928.json`, `gan_600eV_atomic_tube_dense_20260928.json` | At the original VASP/PBE 600-eV, 45.7-GPa contract, the local joint Hessian has a negative direction and sampled 81-point and central 90-point frozen cuts support a bounded mechanism interpretation. The energy–force/stress mismatch and incomplete stationary-point test preclude a strict variable-cell TS certificate. The central cut is not a whole-path or relaxed two-mode surface. |
| Case-specific convergence acceleration | Table 2, `ACCELERATION_EVIDENCE_AUDIT_2026-09-27.md`, Slurm launch audit | Matched-start ABACUS chains use the same first-crossing threshold of 0.10 eV/Å; BTO/HfO₂ counted-launch reductions are 78.5%/38.2%. Negative transfers and path differences remain visible. Do not advertise universal speedup or identical saddle basins. The pre-DFT feasibility gate is a separate safety mechanism; this ablation did not isolate its contribution to launch savings. |
| Installable and usable source package | `docs/RELEASE_NOTES_v0.0.2.md` hosted-artifact smoke; `docs/USER_MANUAL.md`; `tests/test_cli_run.py`, `tests/test_config.py` | PyPI 0.0.2 passes an isolated CLI smoke but predates the later QE/CP2K/ABINIT evidence. The current source is now numbered 0.0.3, which has separately passed isolated local wheel/sdist build, Twine checks, fresh wheel install, and the sdist analytic toy. No 0.0.3 tag or publication exists; neither artifact smoke certifies a user's DFT executable. |

`MANUSCRIPT_EVIDENCE.md` and `FIGURE_LOGIC_AND_STYLE.md` are the current
figure-to-source maps. A pretty interpolated pixel is never an additional
DFT observation. The claim limits above must remain in the abstract, results,
captions, and conclusions when wording changes. The direct manuscript-number
regression `tests/test_cpc_numeric_claims.py` recomputes headline GaN barriers,
BTO blind errors and monotonic endpoint rise, Slurm launch reductions, and
the HfO$_2$ per-cell conversion from committed source tables.

## Technical QA snapshot, 2026-09-29

On the shared working tree, `python -m pytest -q` completed with 616 passed,
one skipped, and 209 warnings (mostly ASE/spglib/Phonopy deprecations).
The clean `17a4f08` Git archive separately passed 597 tests and compiled its
15-page manuscript without missing figures or citations. Increasing the
permitted double-column top-float occupancy then packed GaN Figs. 5--8
beside their discussion across pages 9--11; a fresh working-tree TeX Live
build is 14 A4 pages. Its late pages were rendered and checked for clipping;
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
one skipped). The count is lower than the 616-pass shared-tree snapshot
because the latter includes untracked tests from another work stream. The
archived source test is a reproducibility check, not a fresh DFT audit or a
final CPC/Mendeley deposit.

## Remaining gates for this CPC submission

1. Replace author-order, affiliation, corresponding-author, CRediT,
   funding/computing-acknowledgement, and competing-interest placeholders
   with author-approved statements. Confirm the CPiP article route and
   resolve the Program Summary's pending CPC Library field through the actual
   CPC/Mendeley submission workflow.
2. Obtain a final scientific read-through of every literature comparator,
   pressure and formula-unit normalization, figure caption, and Data
   availability statement. In particular, retain the BTO PBEsol distortion
   versus NEB distinction and the GaN/CP2K and TS limitations. The first
   DOI/source-page pass and five corrected author records are documented in
   `CITATION_AUDIT_2026-09-28.md`; this does not replace the final read-through.
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
