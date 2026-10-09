# E032: paired fixed-control stationary response

This case tests a necessary implementation step of the registered response
comparison, not independent HfO2 evidence, a new harmonic theorem, a TS
certificate or complete B2--B5 freezing. It uses 0 DFT / 0 calculator calls and
does not alter any live production input, code archive or job.

The [API and formulas](../../../../docs/STATIONARY_GAP_RESPONSE.md) retain
reference relaxation offsets and separately scale the initial and bottleneck
charts at the same physical control. Raw centre gap, stationary-anchor gap
and parameter response remain distinct. Frozen-space residuals and clamped
reactions are preserved; caller-declared domains are not automatically proved.

Reproduction uses the saved `autoresearch.sh` with a new output file. Original
local, clean archived and HF reports and delivery receipt are retained after
validation; none is material accuracy or proof of label blinding. See the
dated validation receipt for the actual tested source/archive and results.

## Validation completed on 2026-10-09

- Focused tests: 112 passed in 1.97 s, including 46 new cases. The first run
  had 109 pass / 1 fail: a test guessed the omitted-gradient bound >0.1, but
  the saddle's independently computable value was 0.071. The corrected test
  compares the direct gradient exactly; no physical threshold was changed.
- Clean staged tree `ccfdb0b5` excludes all user-dirty and untracked work.
  The full suite passed: 1162 pass, 2 skip, 313 warnings, 133.93 s; JUnit has
  0 errors/failures and confirms 46 new cases. Source/test bytes are unchanged
  after testing. Later edits are documentation and delivery evidence only.
- Eight analytic groups / 24 points compare separate stationary solves and
  finite gap derivatives. Clean gap and response residuals are <=2.78e-17;
  first/second derivative residuals are <=1.34e-12 / 1.46e-11 respectively.
  These residuals measure implementation agreement, not DFT accuracy.
- The exact archived source was replayed on HF with the existing ICU Python:
  196 float fields agree within 2.09e-12, 85 nonfloat/source-hash fields match
  exactly. NumPy versions remain explicit (clean 1.23.5, HF 1.26.4). No
  software was installed and no remote full-pytest result is claimed.
- All three original reports remain: the working report binds the legacy
  `relaxed_curvature.py` LF bytes, whereas the canonical clean/HF reports bind
  its unchanged Git CRLF bytes. The E031 newline audit is not waived; the
  clean/HF raw-source comparison is strict. That legacy module was not edited.
- At 10:11:33 CST, both live G1 handles remained RUNNING on hfacnormal01:
  preserving 28319570 step 55 / 0.182393 eV/A on node11; reversing 28319571
  step 18 / 0.343628 eV/A on node26. PO->M 28298794 remains COMPLETED 0:0.
  These rounded log observations are not frozen material predictions.
- The production source archive remains SHA256 `df4c12eef03d00d02c3ae357e3f47143b6500b96b7bd12f479f2e2083abd5220`.
  No job/input/optimizer/cutoff/KPT change, extra SCF, early G2 or holdout
  occurred. Actual G1 completion, G2/G3 data, B2--B5 freezing, independent
  advantage and a JCTC-ready paper remain incomplete.

Hashes, source locations, comparison counts and test receipt are in
[`validation_delivery.json`](validation_delivery.json). Retain this necessary
pairing step; do not call it a material result or replace the registered gates.
