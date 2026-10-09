# E044: bounded switching-path Berry measurement

This is the missing G1 property measurement, not a new NEB channel or a
new peak search. Reuse the original preserving step69 and reversing step45
bands (seven interiors each) and three already computed PO endpoint audits.
All nine R3 directions retain positive z with transverse fraction <2e-5.
The preregistered common-axis gate is 1e-3; a failed gate would stop this
one-component experiment, not silently expand it to three axes.

Budget: fourteen output-only SCFs and forty-two fixed-charge NSCFs, one
node/32 genuine MPI ranks per point, maximum two concurrent allocations.
25-minute per-point cap is a resource cap, not a convergence prediction.
SCF preserves the six historical files except for the already validated
`out_chg 1` and `out_bandgap 1` output delta. Original energy/force/stress
reproduction gates are 1e-5 eV/cell, 1e-4 eV/Angstrom and .02 kbar; these
are property reproducibility gates, **not** endpoint relaxation criteria.
NSCF uses the same converged charge, native ABACUS3.10.0/f7cb1d3 Berry
recipe and 2x2x{2,4,8} observable quadrature. No NSCF energy enters a barrier.
100Ry, complete10auDZP, PBE, pseudopotentials and baseline SCF2x2x2 stay fixed.

The new case adapter uses the existing native SCF/NSCF auditors and executor,
without weakening endpoint stationarity/common-cell preparation or the NEB
factory. All source inputs, logs and original E/F/stress are checked first;
cached endpoint SCF/NSCFs are freshly replayed on HF before submission. New
properties live in fresh namespaces; source calculations are never modified.

Variable-cell polarization is tracked in reduced units P/(e|R3|/V), retaining
the explicit native spin-paired period2. Never halve the printed P or choose
the smallest absolute value. The nearest sampled lift starts in an explicit
integer gauge0; half-period ties and quadrature-overlapping alternatives
abstain. Even a unique sampled lift is conditional: nine samples do not
certify unsampled winding, full-BZ insulation, spontaneous polarization,
unmeasured transverse components or a stationary TS. The longitudinal .01
C/m2 gate is a measured sensitivity criterion, not an all-error bound.
No automatic extra sampling, restart or G2 submission follows the summary.

Commands (prepare/run require HF source paths; blueprint is zero-DFT local):

```sh
python -m scripts.hfo2_switching_path_polarization blueprint --root /new/blueprint.json
python -m scripts.hfo2_switching_path_polarization prepare --root /new/batch
python -m pytest -q tests/test_polarization.py tests/test_hfo2_endpoint_polarization.py tests/test_hfo2_switching_path_polarization.py
```

Tests containing synthetic SCFs are protocol checks, not physical evidence.
The following submission record is historical. Completed results and the
delivery receipt below are kept separate rather than revising its meaning.

## Actual submission

Array28418996 was submitted at19:03:04 CST after an empty queue check,
full1285-pass/2-skip regression of the exact immutable runtime archive and
successful HF raw-source/three-endpoint replay. `submission.json` binds the
archive, manifest and scheduler recipe; no live source is overwritten.
At19:25:55 the first eight tasks had completed normally; later tasks remain
bounded by array throttle2. JobArrayTaskLimit is that intentional concurrency
limit, not a failure. Results/whole-path conclusion are not inferred from
submission or partial completion.

An observed pair of concurrent32-rank jobs on one node had disjoint allowed
CPU sets, not overlapping singleton processes. Integer h/min/s native timing
must be parsed as a full tuple (not the leading hour0); it is distinguished
from Slurm allocation seconds and not presented with invented precision.
The finite source/hash/charge audits also incur real launch/I/O cost.

## Completed measurements and interpretation

All fourteen tasks of array 28418996 completed with exit code 0:0. There
were fourteen SCFs and forty-two NSCFs, no new endpoint/peak calculations,
and no repeated DFT after analysis fixes. The output-only SCFs reproduce
the baseline E/F/stress with zero difference at the parsed raw-log precision;
that precision is not a bound on the total physical/numerical error.
The actual six-file contract and SCF charge hashes were checked on HF.

| Candidate | Maximum 224-to-228 change (C/m²) | Minimum half-period margin | Conditional reduced increment |
|---|---:|---:|---:|
| Preserving | 0.000637679742545 | 0.729131419382 | 1.623691304610 |
| Reversing | 0.000838430772605 | 0.611963804869 | 1.623691304610 |

Here 224/228 denote $2\times2\times4$/$2\times2\times8$ longitudinal
quadratures, not different SCF k meshes. Both sensitivity gates pass. Both
nine-image nearest-sample lifts are unique under the declared finite-sample
rule. Their raw-class jumps are modular, not physical jumps. Their equal,
noninteger endpoint increments do **not** establish two winding sectors or
a quantized pump. A unique sampled lift still does not establish unsampled
winding, full-BZ insulation, an absolute spontaneous polarization, transverse
electronic components or a stationary TS. It does not automatically close G1
or authorize G2. The finite nearest-link hypothesis is explicit.

Native integer h/min/s reports sum to 3155 s (28.044444 core-hours at 32
ranks); Slurm allocation times sum to 4029 s (35.813333 core-hours). The
latter includes launch, copying and audit overhead. Neither is substituted
for a full-precision subprocess timer. Concurrent jobs observed on node11
had disjoint 32-CPU allowed sets, not repeated serial processes.

`completed_HF/file_inventory.json` binds 395 exported observable files.
Fourteen raw SCF logs and forty-two raw NSCF logs/band tables replay exactly
to both `offline_replay_HF.json` and `offline_replay_local.json`. Licensed
pseudopotentials/orbitals and bulky charge cubes are not redistributed;
their hashes were checked on HF, not asserted from these omitted exports.

The original production source remains untouched. Its post-DFT summary
encountered a NumPy boolean JSON serialization error, so a new immutable
analysis source was tested and used with zero new DFT calls. An exact known
CRLF/LF native-adapter pair is recognized for this historical archive;
unknown hashes still fail. A separate figure-export scalar serialization
issue was fixed and its incomplete export retained outside Git. The receipt
records these issues rather than disguising them as physical DFT failures.

## Offline reproduction and manuscript figures

From the repository root, choose a fresh output filename/directory:

```sh
python -m benchmarks.hfo2_channels.20261008.switching_path_polarization_20261009.check_export audit --root benchmarks/hfo2_channels/20261008/switching_path_polarization_20261009/completed_HF --output /new/replay.json
python -m scripts.plot_hfo2_path_Berry --output /new/figure
python -m pytest -q tests/test_hfo2_path_Berry.py tests/test_hfo2_switching_path_polarization.py tests/test_hfo2_switching_polarization_export.py tests/test_hfo2_terminal_network_figure.py
```

The [Berry figure](../../../../paper/VARNEB_JCTC/figures/hfo2_path_Berry_20261009/README.md)
contains eighteen plotted rows (fourteen new interiors plus reused endpoints,
with the common initial well repeated), not eighteen new SCFs. The
[terminal-network figure](../../../../paper/VARNEB_JCTC/figures/hfo2_G1_terminal_20261009/README.md)
uses forty existing records and zero DFT calls. Both have source CSVs,
editable exports and actual pixel-review QA. Updated main text preserves the
Pbcn stability/escape-coverage caveat and does not claim G2/G3 or independent
prediction results.
