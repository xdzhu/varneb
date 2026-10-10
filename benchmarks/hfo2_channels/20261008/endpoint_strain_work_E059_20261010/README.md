# E059: measured training well response before interpreting barrier response

This zero-DFT milestone uses the ten existing screened endpoint representations
at the two registered training strains0/+1%. It does not calculate an
intermediate structure, unseen condition, path barrier or additional Hessian.
It implements the already registered B1 controls, not a new physical theory.

| Endpoint representation | Actual0-to1% shift (meV/HfO2) | Imposed-plane trapezoid defect | Full configuration-chord defect |
|---|---:|---:|---:|
| PO+ | +6.806502 | +0.313300 | +0.068647 |
| T descendant | +8.592411 | +0.346594 | +0.204178 |
| M | -15.399977 | +0.018899 | -0.000996 |
| PO-, T pattern preserving | +6.841761 | +0.195524 | +0.090502 |
| PO-, T pattern reversing | +6.841761 | +0.195524 | +0.090502 |

Defects are observed energy-change minus endpoint-work trapezoid, in
meV/HfO2. They are not measured rigorous error bounds. The second diagnostic
includes force and open-cell residual work along the registered endpoint
configuration chord. Intermediate gradients/branches are unmeasured, so
neither trapezoid proves complete derivative consistency or interpolation
accuracy. The+1% T seed is the existing metric-lowered descendant (Ccce),
not a symmetry-restored tetragonal bulk phase.

The imposed-plane work uses the full row-cell matrix, ASE tensile stress
and the explicit biaxial cell direction. It is a local partial, not a
stationary-branch envelope derivative at the nonzero screened residuals.
Analytic tests independently differentiate an energy with mixed atomic/cell
motion and fixed pressure; oblique rotations and shear verify conventions.
`vcneb.strain_work` itself is calculator-independent and calls no calculator.
See[API/definitions](../../../../docs/STRAIN_WORK.md).

## What the registered endpoint-only null controls already imply

At the **seen** upper training condition, all four B1_fixed response increments
are -6.806502meV/HfO2 because the shared initial well rises. If the absolute
bottleneck follows the final well instead, B1_follow gives:

| Registered channel | B1_follow delta B (meV/HfO2) |
|---|---:|
| PO to T | +1.785908 |
| PO to M | -22.206479 |
| Both registered flip correspondences | +0.035259 |

These are assumptions to test against eventual matched barriers, not observed
barrier changes. M stabilization relative to PO+ alone accounts for a
-22.206479meV endpoint difference change; it does not establish a lower
escape saddle or a favorable switching/retention window. A model must
separate actual bottleneck response from these shared-well shifts before
claiming a new mode mechanism or predictive advantage. No channel is selected
here and no+0.5% geometry, label or forecast is created or read.

## Evidence and replay

`endpoint_work.json` preserves raw source/geometry/audit/log hashes, both
endpoint residuals, stresses converted to physical work, ordered-lift source
records, chord components, actual energy changes and all four B1 increments.
The physical INPUT/KPT/pseudo/orbital bytes and original100Ry/full10-auDZP
contract are unchanged. Licensed pseudo/orbital bytes are not redistributed:
offline extraction verifies available raw exports and the earlier HF checks,
not a fresh full-asset byte check. Historical data/receipts remain read-only.

```sh
python -m scripts.analyze_hfo2_endpoint_strain_work --output /new/endpoint_work.json
python -m pytest -q tests/test_strain_work.py tests/test_hfo2_endpoint_strain_work.py
```

The17initial new tests passed locally. Further related, clean-delivery and
HF replay evidence is recorded separately when actually completed. This is
not a full-project or HFpytest claim. No parent or dependency job is changed,
noCI/parameter tuning/monitor/release is added and the overall goal is not
complete. The next material gate remains full matched G2 training paths,
then the existing finite G3 and truly prospective strong-control comparison.

Actual independent source replay passed75related cases in9.13s; the original
HF environment reproduces the analysis at23:15:21with a maximum numeric
difference2.84e-14. All eight used runtime files match the tested archive;
all ten original terminal physical file sets, raw logs and32MPI/EFS were
independently read again on HF. No DFT executable launch was permitted.

Two subsequent fixture checks exposed only legacy archive/checkout newline
expectation errors (each75pass/1fail), both preserved. All three legacy file
line sequences are identical. `source_text_checks.json` records exact raw
executed hashes and separately canonicalLFtext hashes; only legacy text uses
the latter comparison. New scientific code and both paper sources still
require exact raw bytes. No source/input is rewritten or normalized.

The working paper's addedTable4 and work discussion compile to14pages;
updated pages12--14are visually inspected with no clipping or table overlap.
This still is an incomplete material study. At23:16parents remainRUNNING:
flip step8/.363645 andMstep13/.289108, no realfailure; meanSCFs179.91/118.85s
suggest remaining segment caps~03:28/00:53Oct11, not convergence promises.
The previously accepted28692775/76dependency jobs are not modified.

Final corrected fixture delivery passes76cases locally (10.510s) and in the
actual independent archive cwd (7.305s), with no failure/error/skip. Tested
tree88b9214d8185dba6ae92d761132459ebf2be5e73 uses archive SHA-256
1b1b7c35b4a9dee19eb309a64f5985b2ed97a03cf79f653eca13f21d5625065a.
Full-project regression in that same clean cwd then passes1515cases/0failure/
0error with2skips (JUnit383.231s). The skips are the unavailable optional
`ase.io.abacus` reader and non-redistributed licensed VASP POTCAR, not hidden
material convergence claims. Existing deprecation warnings are not failures.
Only reports and user-facing documentation are added after these tests;
the scientific runtime and exact HF-replayed manuscript sources do not change.

The final four-handle scheduler receipt is23:30:26CST: both parentsRUNNING,
both registered childrenPENDING/Dependency. A separate actual `scontrol` read
confirms each parent has12hwalltime ending08:07Oct11, later than the estimated
current segment caps. No submission, cancellation, dependency or runtime
source is changed. This milestone does not complete the overall research goal.

E061later corrects the historical-source fixture's working-copy coupling:
[immutable source bundle and full numeric replay](../replay_delivery_E061_20261011/README.md).
All eight actual executed bytes remain exactly archived, including legacy
CRLF. Future paper prose need not retain the historical snapshot's bytes.
Current analysis now checks all3869numeric entries against the actual HF
report, not just B1 energy increments. This change does not modify any old
receipt, scientific result, physical input or running source/job.
