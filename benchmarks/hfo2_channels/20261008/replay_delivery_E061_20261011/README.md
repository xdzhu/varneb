# E061: historical byte evidence and full numeric material replay

The E059 delivery fixture previously compared **current editable** paper and
library files to the exact source bytes of the historical HF run. It would
reject a later prose-only update despite unchanged scientific inputs/code/
numbers. Exact execution evidence must be preserved, but it is a historical
version contract, not a freeze on every future manuscript revision.

At2026-10-11 00:34:21CST, a source-only HF export rechecked the original
E059 receipt SHA256 and all eight actually executed source-file bytes. It
archived exactly those eight members in`source_bundle.zip` without changing
line endings. Source data147582bytes compress to55195bytes. Bundle SHA256:
`5178718515dfdfae3a03b98b67d82b7c30270f587f1dea1c8de284619b8f24c5`.
The existing E059 receipt/outputs and every production source/job stay intact.
The archive has no INPUT/KPT/PP/orbital, charge/wavefunction or licensed asset.
This export is **not** a new HF numerical replay, native-asset audit or HFpytest.

The fixture now independently verifies all eight archived raw SHA256 values,
the earlier LF-canonical supplement, original receipt identity and archive
contents. It does not require mutable current paper/library text to retain
those old bytes. Current code must instead reproduce **all3869numeric
entries** in the original HF material report, including actual endpoint
energies, volumes, imposed-plane partials, atomic/open-cell chord work,
quadrature defects, residuals and the four B1 increments. Numeric key sets
must agree, values must remain finite, and absolute replay tolerance is1e-9.
The measured maximum difference is2.842170943040401e-14. This computational
comparison tolerance is not a DFT uncertainty or barrier error bar.

Six mutation cases prove that partial, atomic, released-cell, chord, well
and B1 changes are rejected. A re-signed bundle with a corrupted source is
still rejected by the original member digests. A temporary monkeypatch
simulates edits to both current manuscript sources, proves the old contract
would reject them, and verifies that historical evidence survives; it writes
no real paper file. Current scientific runtime/DFT parameters are unchanged.

Local84related cases pass6.82s. Further independent archive verification is
recorded only when actually completed. The full-project1515pass/2skip and
original HF replay in E059 are retained as separate historical evidence,
not repeated or relabeled here.

This is a reproducibility/coverage improvement, not a new material barrier,
G3 surface, stationary-point certificate or prediction. All eight training
handles remain subject to the actual job/convergence gates. The ordinary0.10,
ABACUS100Ry/full10auDZP, finite matrix and unopened+0.5%holdout are unchanged.

```sh
python -m pytest -q tests/test_strain_work.py tests/test_hfo2_endpoint_strain_work.py tests/test_prediction_controls.py tests/test_epitaxial_boundary.py
```

Do not replay the HF export script into its already used single-use namespace.

Actual independent full-archive cwd also passes84cases/0failure/error/skip
in8.146s. Tested treef5209fde9ae8626ecc3aceb8f6b0dc509c86e5cc, archiveSHA256
7c9a427401be4fb599dcf3425332e4760e039f15b61d31cfa32c358023bcbfe8;
all six then-present tracked case/test files match raw archive bytes. Later
changes only append validation/docs. 00:43CSTsacct reconfirms both parents
RUNNING32CPU and all six successorsPENDING; no submission/job mutation here.
The five preexisting user dirty files remain unchanged and un-staged.
