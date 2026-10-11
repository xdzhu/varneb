# E093: compatible terminal audit support for two accepted continuations

This is bounded audit preparation, not another physical experiment or a new
material result. The completed G2 panel is still 4/8; neither G3 selection nor
the +0.005 physical holdout is admitted by this work.

## Why a separate support suite

E091 exposed a real frozen-preparer incompatibility. Its corrected parent
audits remain historical and unchanged. The old E075 exporter also predates
the byte-registered E090 wrapper now executing in jobs 28727368 and 28727382.
Mixing that exporter with the new preparer is insufficient for auditing the
children. This case binds all four current audit support modules together,
while replaying each chain through its **unchanged original E055/E054 runtime**.
No source inside either production runtime is patched.

| Actual child | Parent | Channel | Original Python files | Finite updates |
|---|---|---|---:|---:|
| 28727368 | 28709789 | +1% PO to M | 338 | at most 40 |
| 28727382 | 28723655 | zero-strain pattern-preserving flip | 334 | at most 40 |

Each uses one node/32 MPI ranks, unchanged original six physical input byte
hashes, 100 Ry/full 10 au DZP/PBE, P=0, Gamma 2x2x2, ordinary 0.10 eV/Angstrom,
seven active interiors and fixed endpoints. Nine exact seed caches imply zero
fresh initial SCFs and nominally 7 per completed update. Unexpected additional
calls require inspection; they are not silently treated as successful retries.

## Execution boundaries

The deployed adapter is under the new, audit-only namespace:

```text
/public/home/iai806/abacus/agent-runs/20261011-varneb-clamped-child-audit-E093-r1
```

Run each command in a fresh Python process. `--preflight-only` verifies the
unmodified native accounting identity, submitted ancestry/seed environment,
four actually imported support modules and the original runtime archive.
It accepts RUNNING without relabeling it COMPLETED, creates no terminal
artifact, and does **not** certify all SCFs, final convergence or barriers.

```bash
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  /public/home/iai806/.conda/envs/icu/bin/python -B \
  /public/home/iai806/abacus/agent-runs/20261011-varneb-clamped-child-audit-E093-r1/audit_child_HF.py \
  --job-id 28727368 --preflight-only
```

Only after native successful completion, omit `--preflight-only`. Terminal
mode rejects live/failed/waiting jobs before importing material audit support
or reading summaries. The existing terminal auditor separately checks actual
exit-zero/32-rank accounting, every fresh SCF's complete E/F/stress and six
physical hashes, original runtime bytes, the final complete chain and its
same-boundary ordinary residual. The output is exclusive; an existing audit
or a partial failure is inspected, never overwritten or replayed blindly.

An ordinary-converged chain receives no continuation seed. Only a successful
unconverged 40-update cap receives a latest-nine-cache geometry seed, explicitly
without restored FIRE velocities/history. **Preparing a seed does not authorize
or submit another segment.** Held jobs 28709790/28709791 remain untouched.

Neither driver mode contains submission/release commands, a DFT launch, changes
to physics, a CI request or holdout access. Discrete path peaks are not certified
stationary transition states or a measured barrier-error bound. E092's signed
path-work diagnostic remains training-data evidence, not a forecast advantage.

## Evidence

`tests_local.xml` and `tests_local_v2.xml` retain the actual first and final
local related-test runs. Native compatibility reports, archive replay and
validation metadata are added only after those checks actually run. No native
pytest installation or new full-suite pass is inferred from compatibility.
