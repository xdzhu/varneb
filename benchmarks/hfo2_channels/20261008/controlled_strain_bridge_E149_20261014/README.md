# E149: explicit controlled-strain conditional evaluator

The mode evaluator now accepts an explicit biaxial-control chart, reusing the
existing curvature geometry/stress primitive. Its40-coordinate Hf4O8 chart
includes39atomic/open-cell coordinates and one prescribed substrate parameter.
That final coordinate must not be relaxed. A backend-independent covector check
rejects unsafe conditional-release spaces before any calculator callback.

Actual local related validation:200pass,0failure/error/skip. Checks include
oblique/rotated, tilt/normal and P0/nonzero-P EMT finite energy differences;
both conditional optimizers preserving the control at every callback;
cache/resume/budget identity; malformed raw outputs/geometry; fixed/mixed axes;
and18retained HfO2 zero/+1% native frames in one common controlled chart.
The EMT checks are not DFT uncertainty tests. Native frames are replayed,
not recomputed. No reserved physical+.005point is generated or read.

No material G3 centre/mode/basis/face, full Hessian, prediction or new DFT/job
is selected or executed by this interface milestone. Original input settings,
registered two-bottleneck/two-face budgets and strong-reference comparisons
stay unchanged. Independent snapshot/full-suite evidence is recorded only
after actual execution. The full scientific goal remains active/incomplete.

## Actual independent related verification

Tree6884cf4ca3152e19b3da5be32fc6b90ddf1284b8,
archive SHA5526ab7b335544ccf6be16d50b5d47129356ca4685782911c6f7fb42b9e29b1f.
All6412tracked files exist before tests and all six edited code/test/manual
texts match after source-only newline normalization. Related validation passes
200independently,0failure/error/skip,20.451JUnit seconds; local200 takes14.837s.
Both actual reports are retained. No HF pytest or installation was performed.

At18:09:36CST the additional independent full suite is still RUNNING, not a
passed result. Its preserved exec session is17076 and intended report is
E:/TEMP/varneb-E149-independent-full-20261014.xml. [validation.json](validation.json)
records that observation explicitly; this commit claims related validation,
not a completed full regression. The prior E1482700/2result is not reused as
an E149full-suite result. Subsequent results must be attached after execution.

Only evidence metadata/report copies/plan and ledger follow the tested tree;
the six tested code/test/manual texts remain unchanged. The actual HF read-only
check at18:03:10 finds no hfacnormal01 jobs for iai806;28852667 is COMPLETED0:0,
32CPUs. No new material allocation/DFT or G3/holdout action occurred.

Follow-up during E150: the same independent E149 run actually finishes with
2742passed/2expected skips/0failures/errors,314warnings,1171.095JUnit seconds.
The original report is retained as independent_full_tests.xml,
SHA6235178eefde40751bb65c69ccd7c3c95498c503d7b20c390a3f2e073e766bc4.
The earlier RUNNING record above is historical and not overwritten. This full
suite covers E149, not the subsequently added E150 material driver.
