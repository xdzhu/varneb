# E065: full clean-checkout regression after E063/E064

This is software regression evidence, not a new DFT calculation, model
forecast, optimized HfO2 barrier or JCTC-readiness certificate.

The complete committed tree at
`753c1228ab8ef41578fd37a83fdaa48f3d616f5b` was exported with `git archive`
and unpacked into a fresh directory independent of the user's worktree.
Archive SHA256:
`e5bd596a593b8dc59432007e21bde8ffc962eb455de32eefaa37a90b839d4b52`.
Actual test cwd:
`E:/TEMP/varneb-E065-full-d4dac23eb0a44afbb74de91fe2a41e91`.
Imports of `vcneb` and the new prediction evaluator were independently
confirmed to originate there, not from another installed package.

## Actual complete result

The unfiltered `tests` suite collected1627cases: **1625passed,2skipped,
0failures,0errors**. Pytest reported326.04seconds and314warnings.
JUnit records326.026seconds and timestamp2026-10-11T02:24:22.354234.
`full_regression.junit.xml` retains all individual case outcomes.

The two explicit unverified areas in this checkout are:

- `test_abacus_reader_ignores_optional_eigenvalue_parser`: the installed ASE
  does not provide optional `ase.io.abacus`. This skip is not proof of that
  reader's compatibility, and does not negate separately recorded native
  ABACUS production/EFS audits.
- `test_preparation_keeps_original_electronic_inputs`: licensed VASP POTCAR
  is not redistributed in the clean checkout. This exact-asset preparation
  test is unverified here; no proprietary file was copied to make it pass.

Observed warnings include spglib/phonopy/ASE/fontTools deprecations and two
invalid `\pm` Python string-escape warnings in an existing GaN plotting
script. They were not suppressed or claimed repaired. The console summary
count is recorded, but this case does not claim a complete raw-console file.
The full JUnit is the retained per-test result evidence.

## Reproduce the scope

Tests used Python3.10.9,NumPy1.23.5,ASE3.28.0,pytest7.1.2 on Windows.
OMP/MKL/OpenBLAS thread limits were each1, and unrelated third-party pytest
plugin autoload was disabled. From a fresh archive of the identified commit:

```powershell
$env:OMP_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OPENBLAS_NUM_THREADS='1'
$env:PYTEST_DISABLE_PLUGIN_AUTOLOAD='1'
python -m pytest -q --junitxml=full_regression.junit.xml --durations=10
```

This is not a new full HF regression, every possible ASE/dependency-version
matrix, or a package-installation test. User's five pre-existing tracked
changes and unrelated untracked files were excluded/preserved. No library,
driver, scientific input, live source or job was changed. Later milestone
additions are this report/JUnit/receipts and progress documentation only.

At02:32:11CST, both production handles28661019/28722320 remained RUNNING32CPU,
six successors PENDING. M continuation completedstep5/.129359; flip remained
step17/.260776. These are live unfinished observations, not barrier labels.
The original finite G2/G3/strong-model/independent-holdout requirements remain
incomplete; the full research objective stays active.
