# E162 — local-response evidence integrated into the hafnia manuscript

2026-10-15. Outcome: manuscript progress, not completion of the JCTC goal.
No DFT, scheduler write, production-runtime change or held-out-label access
is performed by this integration.

## Argument and evidence

The added argument is that an ordinary-converged chain or a good displacement
projection does not establish equilibration of coordinates omitted from a
reduced landscape. This is tested with completed native calculations, not
inferred from projection coverage. The main Results retain the original G2
barriers and figures, then add the local curvature, strong-reference and
fixed-coordinate release evidence. The Discussion separates these observations
from the still-unmeasured independent material prediction.

| Manuscript location | Actual source | Claim boundary |
|---|---|---|
| Methods 2.4; Appendix D.4 / A6 | E159: 296 new SCFs, 36 internal dimensions plus one prescribed control, two finite steps | Raw nonstationary curvature, not phonon frequencies, error bounds or TS index |
| Results 3.3; Appendix D.5 / A7 | E159: 28 separate reference/center/step joint-complement records | All four Cmma ranks retained; full T optical control is not declared a failure |
| Results 3.3; Appendix D.6 / A8 | E160: 31 new SCFs and four native cache-only optimizer replays | Two local gradient candidates; identity costs censored, preparation additional |
| Discussion | Above observations under unchanged electronic inputs | Conditional stability, unique branch and held-out prediction remain unmeasured |

The energy changes are per twelve-atom cell, not per formula unit and not
replacement whole-path barriers. The two candidate planes have different
original anchors. Near energies do not establish branch merger. A standard
measured-Hessian quasi-Newton initialization is not claimed as a new algorithm.

Sources: [complete curvature](../G3_complete_curvature_E159_20261014/README.md),
[bounded release](../G3_conditional_canary_E160_20261014/README.md),
[live candidate-sheet protocol](../G3_conditional_faces_E161_20261015/README.md),
[manuscript](../../../../paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.md),
[terminology ledger](../../../../paper/VARNEB_JCTC/TERMINOLOGY_LEDGER.md).
No new external scientific claim or citation is added in this integration.

## Verification and presentation

`test_hfo2_manuscript_local_reduction.py` pins the actual E159/E160 JSON bytes
and reproduces all fifteen added rows, including raw sign counts, every Cmma
registration and actual native call seconds. Mutation tests reject altered
curvature, reference and cost rows. A baseline multiset verifies all eighty
pre-existing table lines, nine figure targets and their bytes, thirteen
external citation targets and unchanged wrapper text. Local manuscript links
must resolve. INPUT and KPT here are byte-exact native snapshots; the potential
and orbital payloads remain at their existing native location, with their
hashes in the original audits.

The four manuscript regression modules pass 37 checks locally and 37 in an
independent archive with verified local imports. This is a manuscript-consistency
check, not a new full-code regression or a repeat native SCF audit. The first
independent run found a Windows CRLF-versus-LF wrapper-text fixture mismatch
(36 passed, one failed). The exact LF text matches HEAD. Only that text check
was made line-ending-portable; physical input and native-data hashes remain
byte-exact. The initial failure report is retained.

The unchanged LuaLaTeX wrapper compiles successfully using the existing TeX
Live installation without new packages or enabling unrestricted shell escape
(the installation retains its default restricted mode). The first 29-page
render split the A8 caption. Its scientific explanation was moved after the
table, leaving a short intact caption. The revised 29-page reading PDF was
rendered completely; all four contacts and detail pages 28-29 were inspected.
There are no overfull/underfull or missing-character log diagnostics. This is
an internal reading draft, not a submission-ready or standalone-link-validated
PDF. Exact build and archive hashes are in `validation.json`.

## Remaining work — unchanged scope

Jobs 28862319 (T) and 28862441 (preserving) continue the existing 25-point
candidate-sheet batches under the original ABACUS/PBE/100 Ry/full 10-au DZP
contract, 32 genuine MPI ranks each on hfacnormal01. No additional allocation
or physics change is introduced here. At 01:19 CST they had processed 7/25
and 5/25 nodes with no reported errors; these are gradient candidates, not
certified stable minima or completed faces.

Next: audit completed native nodes and branches, retain failures as holes,
check released stability before calling a surface conditional-minimized,
then perform the already-registered independent error/prediction comparisons.
The physical +0.005 holdout remains sealed, final representatives and forecast
domain are not frozen, and no extra material, Hessian center or full path is added.
