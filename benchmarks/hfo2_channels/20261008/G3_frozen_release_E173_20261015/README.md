# E173 — eight native frozen-R controls and a matched optimizer-arm launch

The pre-registered eight statics28865656 completed normally on hfacnormal01,
one node/32genuineMPI. Original source383Python files, one conservative attempt
per point and all8native E/F/stress records passed the existing E171 native
admission check. That same check replays the original capped point13baseline
and reparses its12old cached SCFs. All six physical inputs remain unchanged:
PBE/100Ry/full10auDZP/Gamma2x2x2/P0/E0, no climbing image or parameter repair.

`controls_summary.json` and `admission_complete_native_preflight.json` retain
exact original bytes, hashes f49c662d and8f0635ee. The cached original raw
centers and released origins are not recomputed. All4axial Q controls per
center are included; no outcome-based point selection occurs.

## What the ablation distinguishes

For each original-center quadratic and both original probe steps, compare
the frozen-R0 shape and released shape against the corresponding DFT values.
Only the released shapes subtract their already seen released origins.
The affine gradient is retained; no coefficients are fitted.

`released shape error = frozen shape error + missed release-shape change`.

The last term is the model's Q-dependent release change minus the DFT
candidate's Q-dependent release change, each relative to its own origin.
The identity is evaluated point by point, with signed terms. A positive
DFT released-minus-frozen difference is not clipped or silently repaired.
Maxima in separate columns may occur at different points and must not be
subtracted to infer a decomposition. Full pointwise results are in analysis.json.

All energy differences below are meV per twelve-atom Hf4O8cell, not per HfO2.
Divide by4for the latter. These are seen-training discrepancies, not errors
on independent tests, measured physical uncertainty or whole-path barriers.

| Center | Probe step (A) | Max. frozen shape discrepancy | Max. released shape discrepancy | Max. DFT release-shape change |
|---|---:|---:|---:|---:|
| Preserving | 0.01 | 0.075880 | 4.671404 | 4.606959 |
| Preserving | 0.02 | 0.116535 | 4.720758 | 4.606959 |
| T | 0.01 | 0.416264 | 1.327177 | 1.168761 |
| T | 0.02 | 0.403160 | 1.346071 | 1.168761 |

For preserving Q=(-.1,0)A at the first probe step, the signed frozen error
is-.064445, the missed release-shape term is-4.606959 and the released error
is-4.671404meV/cell. At T Q=(+.1,0)A the signed terms-.416264and+.378430
partly cancel, leaving-.037833. Neither favorable cancellation nor the
small frozen errors are selected as independent model success.

These matched observations locate most of the preserving four-point shape
discrepancy in the release-response term. They do not isolate anharmonicity,
finite-step effects, finite released residuals or branch changes as a unique
microscopic cause. Positive raw model H_RR does not certify the relaxed
DFT branch. R=0 does not freeze the entire cell because Q contains cell motion.
This ablation is not the full B2-B5/T/Cmma strong-control prediction matrix.

The statics cost8newSCFs/8attempts,8.136814completed-call proxy core-hours;
the complete scheduler allocation is1021s*32/3600=9.075556core-hours.
The earlier296curvature SCFs, conditional candidates and whole-path
preparation remain additional. No net acceleration factor follows.

## Bounded matched arm

After the8controls complete and pass native audit, the unchanged E171
one-shot driver submits28866132, validates its held allocation and releases
it. It uses1node/32MPI/1h, at most12new attempts, and the already reserved
E170point13start/Q/seed/.5A total domain/.03released-gradient target/.05trust.
Old exact cache hits are separate from new SCFs. The old12SCF failure is
censored, not a denominator for speedup. No old runtime is overwritten.

The submission/native receipts are archived here. They prove admission and
input identity, not material convergence or a beneficial optimization result.
Pending/running counts remain part of the two-allocation cap. No additional
training reservation is charged: preserving30/81,51left;T29/81,52left;
all8independent slots per face and the physical+.005holdout remain untouched.

## First new native SCF and a provenance-path fix

At05:11CST the new arm is RUNNING with4complete new SCFs/5reserved attempts,
no observed backend failure and no final arm summary. The first new structure
exactly matches E170's zero-DFT novel proposal. Independent native E/F/stress,
full gradient, fixedQ/control, six assets and32MPI checks pass; its released
gradient norm.069855eV/A is not convergence. This is not a speedup result.

The first read-only audit initially rejected its validation dictionary even
though energy/gradient/stress/Q/domain all matched. The sole difference was
the spelling `bundle/../arm` versus `arm` for the same raw directory. A second
diagnostic also exposed a NumPy-boolean JSON-serialization error. The auditor
now first resolves both paths strictly to the exact expected native directory
inside the permitted arm, then revalidates through the recorded path spelling.
Every validation field still compares exactly; another real directory or an
escape from the allowed root is rejected. A regression test covers this case.

No physical tolerance, raw result, source worker or allocation was changed,
and the job was never stopped/restarted for this metadata-only mismatch.
`first_native.json` preserves both path spellings and the successful raw
validation. The raw geometry radius.49999999999999994A already passed the
unchanged original guard; no numerical-domain tolerance was loosened.
