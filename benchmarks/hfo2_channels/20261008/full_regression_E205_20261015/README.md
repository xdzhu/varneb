# E205: complete offline regression after the approved manuscript reorganization

This checks the existing software and evidence tests while the separately
submitted released-curvature job runs. It adds no DFT calculations, test
coordinates, material cases, physical settings or prediction claims.

## Initial complete run and the actual defect

A clean Git export of `4d03befe07fc053b2012097337c6dd180517ca99`
collected 3506 tests: **3502 passed, 2 failed, 2 skipped**, with 314 warnings.
The original XML is retained as `first_full_tests.xml`, not replaced by the
successful targeted check.

Both failures came from `tests/test_hfo2_E094_manuscript_work.py`, which still
located the old Results 3.4/3.5 after the author-approved four-part rewrite.
Because the old ending heading no longer existed, its section helper read
through the rest of the manuscript. This could make some row/prose checks
pass using material outside the intended section. The Cmma range check
actually consumed 139 numbers rather than the intended 16; the old section
number assertion also failed. This was not a changed DFT result.

## Bounded repair

The same work evidence now resides in Appendix E.1, ending at E.2. Both
headings must occur exactly once and no intervening heading is accepted.
The eight source-derived signed work rows, four prose partitions, all four
registered rank-four Cmma frames and all 16 range endpoints retain their
original numerical tolerance (`5.0001e-7`, `rtol=0`). Section/link checks now
cover the approved four Results and the two actual appendix locations.

Six added tests reject missing/duplicate boundaries, borrowing a missing row
from elsewhere, and each of the 16 altered Cmma bounds even when correct
values remain outside the target section. The targeted module passes all
21 tests. Only this test module changes executable source; production code,
the manuscript, native runtimes and data remain unchanged.

## Full verification

The second clean export contains the exact revised test module and otherwise
the same committed source. Its complete regression collected 3512 tests:
**3510 passed, 0 failed, 0 errors, 2 skipped**, with 314 warnings. All 3506
original test identifiers are retained; the only additions are the six
boundary/mutation tests above. The two skipped identifiers and reasons are
identical between runs. `final_full_tests.xml` is the full rerun, while
`target_tests.xml` records the separate 21-test working-tree check.

Both source snapshots and all temporary fixtures are on C:, not the nearly
full E:. Runs use the existing Python 3.10.9 / pytest 7.1.2 environment, one
numerical thread, existing local data and simulated calculators. No package
installation, scheduler submission, CI service or licensed DFT execution is
part of the test command:

```text
python -m pytest -q --tb=short --disable-warnings --durations=15 \
  --basetemp=<fresh-C-scratch> --junitxml=<fresh-report.xml> \
  -o cache_dir=<fresh-C-cache>
```

The original two skips are explicit coverage limitations: the environment
does not provide optional `ase.io.abacus`, and a clean public checkout does
not redistribute the licensed VASP POTCAR needed by one preparation check.
They are not counted as passes and no dependency is installed to hide them.
Warnings are not a claim of warning-free execution.

Source trees, archive hashes, XML counts and protected-file hashes are
recorded in `validation.json`. The revised executable-source tree is
`1777a194595cb0794852c5431aafbe0ac4e4ae5c`; its Git-normalized revised test
blob matches the working-tree test exactly. This passing software regression
does not establish released-branch stability, a TS,
independent prediction accuracy, scientific novelty or completion of the
JCTC study. Those material-data gates remain separate.

At 12:28:26 CST, the existing curvature job 28868827 was still RUNNING on
node9 with 32 CPUs: 82/272 completed result records, no failure records.
The mean completed SCF cost gave an approximately18:41 completion scenario,
not a deadline or guarantee. No job was added or restarted in this milestone.
