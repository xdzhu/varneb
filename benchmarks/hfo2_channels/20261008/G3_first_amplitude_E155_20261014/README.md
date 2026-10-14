# E155: first complete amplitude, unstable one-coordinate complements

This is new measured local-curvature evidence at the two existing E150
centers, not a new DFT batch or a certified stationary index. The original
28855614/28855615 allocations continue unchanged. The physical+0.005holdout,
conditional-plane choice, forecast freeze and ordinary NEB threshold remain
untouched. No source/SCF/input repair, ASR operation or new author-mode array.

## Native evidence and independent arithmetic

The unchanged E152 auditor was executed in its separate native namespace
against the immutable E150 runtime. It reparsed an ordered prefix of76T and
77preserving probe SCFs, original cached centers, all six actual input hashes,
genuine DSIZE32, electronic convergence, ordered geometry, raw E/F/stress and
the full physical gradient. No failure artifacts were present. These153
replayed probe records include earlier audited records; they are not153new
DFT calls triggered by the audit. Conservative reservations and SCF-directory
counts can exceed the snapshotted complete prefix while the jobs are running.

Both prefixes contain all37signed pairs at.01Angstrom:36translation-free
internal directions and one separate prescribed-strain control. Only one pair
at.02Angstrom is complete in each prefix. Neither second matrix is complete.
Native report bytes were copied unchanged, not normalized or regenerated:

| Report | SHA256 |
|---|---|
| T_native.json | 1cc4269a8f1437acddb94c9a7c264f52bae888669db2bcc36521982d280666ea |
| preserving_native.json | 9bfeb3c8a7eca7533246ffd61e6e43df17262fb8ef274dfa87965c1301a10efa |

`scripts/summarize_hfo2_G3_curvature.py` is an offline arithmetic replay of
these hash-pinned native audits and the unchanged E150 manifest. It does not
claim another raw-SCF reparse. Only complete blocks receive a spectrum; the
external control is excluded from every internal instability count. The
symmetrized internal matrix and its restriction perpendicular to the native
path's measured internal tangent are diagonalized without a physical sign
threshold or a pseudoinverse. All tolerance checks are algebraic, not a
measurement of physical DFT uncertainty.

## Actual first-amplitude observations

Curvatures are eV/Angstrom² in the declared joint metric, not frequencies.
Gradients are eV/Angstrom; a full internal norm is not the ordinary NEB
maximum projected image force and does not replace its.10criterion.

| Quantity | T-derived escape center | Preserving switching center |
|---|---:|---:|
| First internal eigenvalue | -0.276193429 | -1.271818849 |
| Second internal eigenvalue | -0.164648089 | -0.447449785 |
| Third internal eigenvalue | 0.053215645 | 0.173653492 |
| Lowest tangent-complement eigenvalue | -0.176131009 | -0.447450308 |
| Raw internal / tangent-complement negative counts | 2 / 1 | 2 / 1 |
| Measured internal gradient norm | 0.108189564 | 0.305934312 |
| Prescribed-control gradient | 0.670977160 | -0.020117922 |
| Normalized internal tangent / lowest-direction absolute overlap | 0.832681934 | 0.998560997 |
| Raw reciprocity relative defect | 0.000406464 | 0.000272021 |
| Symmetrization operator change | 0.015143022 | 0.009867690 |
| Maximum energy/gradient diagonal disagreement | 0.108566099 | 0.189041460 |

The control gradient is a substrate reaction derivative, not a failed
relaxation requirement. Neither physical internal center is stationary.
At a nonstationary geometry, Hessian signs depend on the declared chart;
these are not second-order saddles or resolution-certified instability counts.

The native unit path tangent has excluded translation norms0.101688097T and
0.007027557preserving. The arithmetic diagnostic explicitly reports this
gauge contribution, projects the tangent into the measured internal space
and normalizes there. The source path and raw forces are never changed.
Consequently the normalized overlap above is deliberately not the E152
unnormalized native-tangent contraction. The initial helper incorrectly
required the native tangent to have no translation component; it refused
before writing a report. That analysis guard, not any material input, was
corrected to expose the projection and excluded norm explicitly.

## Consequence, not a new claimed method

If the measured symmetric matrix has two negative directions, no choice of
only one retained internal coordinate can leave a positive eliminated
complement: its codimension-one complement must intersect the two-dimensional
negative subspace. The actual tangent complements also explicitly retain one
negative direction. This elementary dimensional fact is not a novel algorithm.
It warns against treating a single internal mode plus external strain as an
already stable conditional surface and blindly applying stable-block Schur
elimination or unconstrained complement minimization.

Retaining both negative eigenvectors would leave a34-dimensional positive
spectral complement **for this one measured matrix**. That is only an
algebraic diagnostic, not the final chosen axes, branch or a stability
certificate. This complement still has gradient norms0.107972638T and
0.305934117preserving, so the centers are not conditional minima either.
The second amplitude, numerical consistency, finite-domain branch tests,
strong fixed-reference comparisons and independent forecast tests remain.
Do not silently add a third center, a3Dsampling budget, new thresholds or new
electronic parameters in response. Any final plane/branch choice remains
within the existing finite protocol and must be made before holdout labels.

## Reproduction and scope

```text
python -m scripts.summarize_hfo2_G3_curvature
  --manifest benchmarks/hfo2_channels/20261008/G3_curvature_E150_20261014/manifest.json
  --T-report benchmarks/hfo2_channels/20261008/G3_first_amplitude_E155_20261014/T_native.json
  --T-sha256 1cc4269a8f1437acddb94c9a7c264f52bae888669db2bcc36521982d280666ea
  --preserving-report benchmarks/hfo2_channels/20261008/G3_first_amplitude_E155_20261014/preserving_native.json
  --preserving-sha256 9bfeb3c8a7eca7533246ffd61e6e43df17262fb8ef274dfa87965c1301a10efa
  --output FRESH_EXTERNAL_SUMMARY.json
```

The output must be fresh; pin mismatch and overwrite are refused. Existing
manifest/source/native bytes remain unchanged. Tests cover the actual case,
internal orthogonal reparameterization, external-control exclusion, explicit
translation gauge and malformed/partial/failure reports. Independent-checkout
results and exact scalar CLI replay are recorded in validation.json after
execution. No full-suite, H1/H2, TS or manuscript-readiness claim follows.
