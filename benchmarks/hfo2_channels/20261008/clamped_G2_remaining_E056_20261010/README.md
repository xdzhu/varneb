# E056: the remaining six registered clamped training starts

This is a zero-DFT preparation and delivery milestone, not six new converged
paths. Together with the two E053 zero-strain pilots, these starts complete
the **eight-start** finite G2 matrix. Two existing continuations remain live;
none of the six starts in this directory has been submitted.

## Actual preparation

| Training substrate | Remaining registered channel | State |
|---|---|---|
| 0 | PO to T | HF geometry and endpoint-cache preflight passed |
| 0 | PO flip, T pattern reversing | HF geometry and endpoint-cache preflight passed |
| +1% | PO to T | HF geometry and endpoint-cache preflight passed |
| +1% | PO to M | HF geometry and endpoint-cache preflight passed |
| +1% | PO flip, T pattern preserving | HF geometry and endpoint-cache preflight passed |
| +1% | PO flip, T pattern reversing | HF geometry and endpoint-cache preflight passed |

HF completed the preparation at 2026-10-10 21:18:05 CST, using the already
tested immutable source commit `86d2637ebcd6c502df43b91ca03977fd50120716`.
All 338 Python files in its vcneb/scripts/examples trees were compared
byte-for-byte against the original tested archive before execution. The
archive SHA256 is
`cbe6552a3fa4237fec2098c08c8056735a9a95df4e18691d09fed6b4f24a613c`.
Its earlier clean full regression passed 1459 tests, with 2 skipped.

Each chain has nine total images: seven uncached moving images and two
fixed, exactly matched endpoint caches. All 54 geometries passed periodic
lift, ordered-atom and clamped-boundary preflight. All 12 cache accesses
checked the original six physical file bytes, raw-log pins, native 32-rank
execution and full energy/forces/stress. The ten previously screened
endpoint representations were reused; no endpoint relaxation was rerun.
G1 supplies only its registered ordered geometrical mechanism, never
free-cell energy, force or stress labels for the clamped ensemble.

The substrate fixes cell rows 0 and 1; the third vector and atoms remain
movable, at P=0. All retain original ABACUS 100 Ry, full 10-au DZP,
original INPUT/KPT/PP/orbitals, 32 true MPI ranks and thread count 1.
Ordinary NEB remains 0.10 eV/Angstrom, without CI. Original physical
hashes and exact raw-source identities are retained in each manifest
and endpoint-cache audit. Licensed PP/orbital bytes are not redistributed.

## Contents and evidence

- `prepare_and_preflight_HF.sh`: the executed preparation, guarded against
  reusing an existing output namespace. It starts neither a Slurm job nor DFT.
- `seeds/`: six calculator-free trajectories, endpoint POSCARs, manifests,
  exact fixed-endpoint cache records and HF geometry preflight receipts.
- `geometry_preflight/`: validate-only reports and initial-chain snapshots.
- `cache_preflight/`: the actual 12 endpoint-cache audit records.
- `executed_source/`: exact historical preparation/helper bytes, which may
  differ in line endings from a subsequent Windows checkout.
- `preparation_receipt.json`: original HF receipt, including 338 source
  checks, per-chain manifests and raw cache hashes.
- `local_tests.xml`: 29 tests passed, zero failures/errors/skips, in
  104.664 seconds: the eight new material-delivery checks plus the 21
  existing clamped-chain tests. This is local pytest, not HF pytest.
- `worktree_related_tests.xml`: all six selected material/response/cost
  test files pass, 94 cases in 107.11 seconds. The first command remained
  in the working checkout despite creating an archive; it is not called
  a clean-archive test. A separate test with explicit clean-directory cwd
  is required for the delivery check. This harness issue launches no DFT.
- `clean_related_tests.xml`: the separate actual archive-cwd replay passes
  all94cases in107.75seconds, with zero errors/failures/skips. The tested
  staged tree is `2487220b4a696bbdcdf0022588229f9e793a5141`; all76case files
  present at that test, and both manuscript text sources, were independently
  byte-compared to the archive. Later receipt additions change no executable,
  test, scientific source, geometry or physical input.
- `snapshot_live_jobs.sh` and `live_jobs_snapshot.json`: one fresh actual
  scheduler/log observation, not a new monitor or raw material TS audit.
- `validation.json`: bounded delivery scope and historical bundle pins.

The exported HF results tar SHA256 is
`3c6effb9d54f2bfda3affedad287f60867e73602bedef7b0876bd5e98e721cfb`;
the downloaded archive was checked before extraction. No DFT executable
launch was permitted during preparation. The existing HF environment was
not installed or upgraded. The source and running production directories
were not replaced.

From the repository root, replay the material checks without HF or licensed
assets:

```sh
python -m pytest -q tests/test_hfo2_clamped_chains.py tests/test_hfo2_G2_training_matrix_delivery.py
```

## Deferred calculation, not a new monitor

Reconcile `squeue` and, for missing handles, `sacct` before each actual
submission. The study may hold at most two simultaneous hfacnormal01
allocations; unrelated jobs must remain untouched. First finish the two
zero-strain continuations 28661019/28661020, then fill the remaining two
zero-strain channels before the four +1% channels. A continuation remains
the same independent chain, not an extra trial. The archived pilot script
is a ten-step/four-hour segment template, not an automatic convergence
controller or permission to duplicate an active chain.

At21:38:01CST both continuations were RUNNING at32actualCPUs, on node6/25,
with no actual `band/vcneb_failure.json`. Flip step4/log force0.863395;
M step6/log force0.455959eV/Angstrom. Fresh internal SCFs average163.825/
114.497seconds respectively. With seven moving images per update, the
remaining16/14steps imply about5.10/3.12hours, i.e. Oct11around02:45/00:45
to the segment cap. This is not a convergence promise, full raw-chain audit
or budget exemption. No unchanged-state polling loop is added.

The +0.5% held-out condition was neither generated nor read. G3 stays
within the previously registered two-bottleneck sampling budget. Endpoint
screening and starting-path preparation do not establish barrier ranking,
TS certification, a response forecast, material prediction accuracy or
JCTC readiness. Actual matching G2 paths and errors are still required.
No new material, CI, release, PyPI publication or automation is added.
