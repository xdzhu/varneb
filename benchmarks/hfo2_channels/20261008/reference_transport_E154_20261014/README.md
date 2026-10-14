# E154: strong-reference tangent transport, not a predictive gain

The preceding E153 turn made manuscript/numerical replay progress. This
milestone connects the already registered mass-weighted original-T reference
spaces to the two actual E150 local length-valued curvature charts. It does
not select conditional axes, infer new phonons, supply missing curvatures,
freeze predictions, launch DFT or touch the running E150 runtime.

Public API: [atomic reference transport](../../../../vcneb/reference_transport.py),
documented in [the advanced recipe](../../../../docs/ATOMIC_REFERENCE_TRANSPORT.md).
Original source arrays/forces remain unchanged. All returned cell/control
components are exactly zero. The analysis translation restriction and its
discarded-work term are explicit; this is not source ASR or an assumption
that the raw net force vanishes. Rank loss is rejected, never shortened.

## Actual source use, not private-array redistribution

The verifier checks the original E023 registration/author hashes, T/Gamma
bytes and E150 manifest. Seven spaces are retained at each frozen center:
T geometric rank3, two complete optical doublets rank4, all optical rank33,
and separate Cmma frames103/106/114/127 of rank4 each. The stronger complete
T optical bank uses the existing matrix, not new DFT. Different ranks are
reported, not falsely called same-dimension controls. All spaces lie within
the actual E150 declared measured span; that does not mean its SCFs finished.

Only scalar diagnostics and source hashes are exported. Raw and transported
author vectors remain private inputs outside Git. Existing basis/matrix and
endpoint preparation is not made cost-free by this zero-new-DFT conversion.
Frozen center force records are used, not claimed as newly raw-SCF audited.

The first verifier over-required the entire substrate-reference cell to
equal T and correctly refused before output. Only two rows are prescribed;
the third is open. The helper is corrected to check those two rows, with
the original complete manifest still hash-pinned. No geometry, boundary,
force, physical parameter or running source was changed. The generic
implementation also explicitly normalizes the chart's non-unit translation
directions and tests pure-translation rank loss; source modes are untouched.

Actual14source/center pairs preserve all ranks. Maximum reconstruction defect
2.160e-16(1/sqrt(amu)); independent original-T/local atomic work-covector
defect6.939e-18eV/(sqrt(amu)Angstrom); internal+translation identity
2.819e-17 in the same units; measured-span projection defect2.341e-15.
These are algebraic checks, not DFT errors, mode-energy partitions or
physical stability/TS/prediction certificates. No physical+.005structure
or label is generated. The full research goal remains incomplete.

## Reproduction

```text
python -m scripts.verify_hfo2_reference_transport --source-root PRIVATE_ORIGINAL_AUTHOR_DIRECTORY --output FRESH_SCALAR_REPORT.json
python -m pytest -q tests/test_reference_transport.py tests/test_biaxial_mode_evaluator.py
```

Actual local first implementation checks60passed. Further focused and
independent checks, exact source identities, native-job observations and
delivery outcome are recorded only after execution in validation.json.
The preliminary valid scalar export remains outside Git under E:/TEMP;
the final verifier contracts atomic forces directly, never invents a
zero physical cell/control gradient. Arrays and actual running curvatures
are not copied or rewritten by this recipe.

## Actual independent delivery checks

Final focused checks:77local/77independent,0failure/error/skip,
3.905/5.966JUnit seconds. Fresh tree47a1923b9e6c918fe04c97b1d4f4dc6e16383d88,
archive094f9ab21d03a54157c6cc46a09c911575ad32e33a78de73c4ed6e28a5b0066e;
all6481tracked files present before testing and both import roots verified.
Six implementation/test/manual texts remain normalized-equal to the tested
snapshot. The actual private-input scalar CLI replay also runs in that
checkout:all14rows/scalar/input-source fields match exactly. Only two
executed source SHA fields differ by proved code LF/CRLF export. Physical
data and private source arrays are not normalized or redistributed.

The final helper adds a refusal to write audits inside the private author
source tree. Compared with its preceding valid export, only the helper SHA
changes; all scalar/source fields are exact. Both preliminary exports remain
outside Git. A physical source with zero excluded-translation norm is valid;
the test does not require an artificial positive net-translation signal.

20:43:41CST:both original curvature allocations remain RUNNING,32CPU each,
70/148and69/148complete,node10/node27. This status reads receipts only,
not raw EFS again. The earlier successful E153 GitHub push is now separately
confirmed by a successful native remote-main read; no proxy setting changed.
No new job/SCF, live source, physical contract, reference registration,
conditional selection or held-out input. Full research goal still active;
this adapter is not a full curvature, prediction or JCTC-completion claim.
