# Terminal switching health segments: audit and ordinary continuation

The first 10-step segments 28288045/28288063 completed normally, but their
ordinary residuals are 0.980656799/0.700278856 eV/A, not below 0.10.
Each complete step10 contains nine raw, ordered Hf4O8 SCF evaluations. All
18 were matched to snapshot geometry and the byte-fixed physical contract,
then replayed with full energy, force and stress. No SCF was added for export.

The restart folders preserve those geometries and name all exact raw caches.
The actual immutable R12 production factory checked every cached value with
its DFT entry point disabled. The two preflight receipts retain source,
verifier, manifest and trajectory hashes; this is not a mock calculator check.
The actual R12 material CLI also passed geometry/continuous-lift validation.
Prepared geometry continuation uses a fresh FIRE state, not restored momentum.

Jobs **28319570/28319571** were submitted exactly once in a new R17 namespace,
with `afterany:28298794:28300425` to preserve the two-active-chain limit.
Each requests one node/32 MPI/OMP1 on hfacnormal01, at most80 steps/24h,
nine total/seven active images and two cached endpoints. This extends the
same two registered G1 candidates, not the material/strain/channel matrix.
No cutoff, orbital, pseudo, k-point or SCF parameter changes and no CI.
Details and actual Slurm verification are in `../submission_handles_r17.json`.

Reproduce the production-cache audit on HF in a **new** output namespace:

```bash
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
python verify_production_cache.py --production-source /path/to/immutable/source \
  --seed /path/to/audited/restart --output /path/to/a/new/cache-check
```

Raw SCFs are required for this command and remain on HF. Local evidence replay
needs no remote DFT: run `pytest -q tests/test_hfo2_observation_restart.py`.
Git attributes preserve the hashed bytes, including the executed verifier.

## Updated representation evidence, not optimized barrier ranking

`../channel_network/terminal_switch_step10_analysis.json` compares the four
terminal first-segment observations from the same physical PO+ initial state,
preserving the older step1 analysis. All four residuals fail the ordinary
criterion; sampling/error bounds are missing and selectivity remains unset.

The preserving flip's discrete maximum fell from145.779 to97.201 meV/f.u.
and the reversing candidate's from435.169 to410.744 meV/f.u. These are
**unconverged discrete maxima**, not final switching barriers. The preserving
candidate and PO→M observation (95.199 meV/f.u.) differ by only2.002 meV/f.u.;
the current data do not establish a reliable ordering or a favorable window.

Both flip energy maxima remain dominated by the cell block. The preserving
peak has atomic/cell residuals0.397045/0.980657 eV/A; the reversing peak has
0.143499/0.691671. However, the reversing chain's largest residual is now
atomic at another image. Optimization bottlenecks can move between images
and blocks; one early mechanism label must not be imposed on the whole run.

At the peaks the original T-pattern triplet captures57.39% and essentially
0% of parent-relative squared displacement. This reinforces the need to
select actual joint bottleneck directions; it is not an energy decomposition,
local Hessian, TS certificate or proof that the candidates remain distinct.
No cell freezing, model tuning or extra DFT was used to produce these results.

## Final delivery validation

The complete clean archive f575cd8f / SHA2562c00c537 passed867 tests,
2 skipped in84.72s, excluding unrelated working-tree files. It contains all
four observation/restart trajectories; the generic `*.traj` ignore rule was
detected during packaging and those exact artifacts were explicitly added.
The historical seven-input analysis is tested against its original named
inputs, not rewritten to include new observations. Three new tests cover
the terminal flip evidence, exact caches and refusal of a premature ranking.
The executed verifier hash is identical locally, on HF and in the archive.
See `validation_delivery.json`; following changes are receipts/docs only.
