# E055: nested stationary-response controls and honest cost accounting

This milestone supplies backend-independent **local control algebra**, not
HfO2 forecast accuracy or a completed JCTC result. No DFT or held-out label is
generated; no physical settings or running production sources change.
The two original G2 continuations remain separate live calculations.

## What is implemented

- One measured-column dataset constructs frozen, atomic-only, joint release
  and training-instability promotion controls. Unknown complements remain
  unknown; unstable release is not pseudoinverted or assigned a zero barrier.
- Pair a well and bottleneck at the same physical parameter, retaining both
  responses, nonstationary-anchor corrections, signed gaps and local domains.
- Count the entire visible preparation dataset, including failed/unused
  measurements and exact shared-source reuse. Unknown timing/resources stay
  unknown, not free. No prediction error or CPU utilization is invented.

The [technical contract](../../../../docs/NESTED_RESPONSE_CONTROLS.md)
explains the standard harmonic algebra, coverage and failure semantics.
B5 is only the **local increment-dimension arm**, not an anharmonic branch
search. When joint release is stable B4=B5; that equality cannot be claimed
as improved prediction. The registered strong T/Cmma/gauge controls, complete
same-boundary training, branch/probe validation, prospective freeze and
unseen error tests still determine whether the research claim survives.

## Evidence and reproduction

The tested immutable source is86d2637ebcd6c502df43b91ca03977fd50120716;
clean archive SHA256cbe6552a3fa4237fec2098c08c8056735a9a95df4e18691d09fed6b4f24a613c.
The full clean local suite passed1459tests,2skipped,0fail/errors in318.96s.
The first175related checks also passed. JUnit and validation receipts are
archived separately from the runtime source; the user's five tracked edits
and unrelated untracked files are not included in that runtime/test tree.

`analytic_local.json` and `analytic_hf.json` test48paired model points,
rotated charts, distinct control scales, incomplete measurements and one
unstable cell-direction promotion. Independent full-block stationary solves
and finite derivatives agree. Analytic response errors are **implementation
residuals**, not DFT uncertainty or HfO2 prediction errors.

The first HF verification passed the analytic checker then stopped because
the existing ICU environment has no pytest. That failure is retained; no
package was installed, no source/output overwritten and no DFT launched.
The separate follow-up performs the numeric/cost-semantic checks without
pytest and exports only old, completed pilot logs. The full suite remains a
local clean-source result, **not** a claimed HF pytest pass.

An initial local delivery replay passed5/failed1: the unprotected historical
`audit_hfo2_static_replica.py` had LF in the working tree but Git's Windows
archive emitted CRLF. The runtime-source hash therefore differed from the
later working-tree file despite unchanged content. Preserve the actually
executed source bytes under `executed_source/` and compare their exact hashes,
not a different checkout's text. No physical input/hash gate was relaxed and
no SCF rerun. The archived subset fingerprints the measured operations, not
a claim that it is a complete standalone package.

The material-gate helper's first direct-file invocation lacked repository
PYTHONPATH and failed before reading inputs; the explicit repository import
path then passed. A diagnostic `git diff` also used the non-Git clean archive
as cwd and was repeated in the checkout. Both caused zero DFT/file mutation.
The first progress observation checked the wrong optional failure filename;
the separate corrected receipt checks actual `vcneb_failure.json` and retains
the old observation. These delivery failures are not electronic failures.

The HF cost receipt reads all154fresh internal SCFs from28574708/28574709,
checks all original six physical bytes, every pinned raw log, full native
E/F/stress and oneDSIZE32. Actual `sacct` allocation32CPUs is recorded, not
inferred from the requested mpirun string. Earlier endpoint/preparation and
ongoing continuation costs are **excluded**, so this is not total G2/study or
future predictor cost. Transport wall seconds times CPUs is a resource proxy.
Scheduler allocation and SCF cost are different recorded quantities.

The material negative gate sends the actual unfinished clamped step10
observations to the complete-network entry point and confirms refusal before
any prediction namespace is created. Both supplied sources are zero-strain
pilot observations, not claimed0/+1% network labels. This is a typed-input and
coverage refusal, not an accuracy score of an attempted material predictor.

From a repository root:

```sh
python -m scripts.check_nested_response --output /new/analytic-check.json
python -m pytest -q tests/test_nested_response.py tests/test_response_cost.py tests/test_nested_response_benchmark.py
```

The bounded HF exporter is intentionally registered to E053 only:

```sh
python -m scripts.export_hfo2_scf_costs --run-root /public/home/iai806/abacus/agent-runs/20261010-varneb-clamped-G2-E053-r1 --output /new/prior-cost.json
```

It cannot process a live/unregistered run, silently omit a failed call or
launch ABACUS. Source/raw hashes and per-call identities remain in the receipt;
licensed pseudopotential/orbital bytes are not redistributed here.

## What this does not complete

G2 still needs complete, audited competing paths at0/+1% under the common
substrate. Current numerical interfaces do not supply material Hessian
columns, independent curvature error or branch/TS certification. G3's two
bottleneck/two-surface budget is unchanged. No+.005held-out calculation,
new case, CI, release, PyPI publication or new monitor is created.
The goal remains active; prospective accuracy and paper readiness are pending.
