# E156: strong-reference complements on the same measured local matrix

The two E150 batches are still running. This offline study uses precisely
the first-amplitude native prefixes archived in E155, not a completed
two-amplitude Hessian, additional SCFs or held-out data. All seven registered
references see the **same** symmetrized36-dimensional fixed-strain internal
matrix at each center. The external prescribed control is excluded.

## Actual comparison

Source directions are transported with the existing audited tangent/mass/
translation-gauge map. No new parent Gamma calculation, fitted rotation,
rank truncation, source ASR or frequency transfer from LDA to PBE. Keep the
T geometric triplet, both complete T optical doublets, complete33T optical
directions, and all four separate Cmma registrations103/106/114/127.
The four Cmma cases are not merged into16directions or selected by outcome.

For each retained atomic span Q, explicitly restrict the measured H to:

- the complement of Q in the33-dimensional translation-free atomic space;
- the complement of Q in the full36-dimensional atomic/open-cell space.

The table gives the lowest **joint-release** curvature(eV/Angstrom²) and
raw number of negative joint-release directions. Zero is only the raw sign
counting point, not a physical stability floor or error bound.

| Retained reference | Rank | T minimum | T negatives | Preserving minimum | Preserving negatives |
|---|---:|---:|---:|---:|---:|
| T geometric triplet | 3 | -0.230562155 | 2 | -1.271818328 | 2 |
| T two complete optical doublets | 4 | -0.175388281 | 1 | -0.246012653 | 1 |
| T complete optical space | 33 | 2.956213824 | 0 | 3.090797687 | 0 |
| Cmma103 | 4 | -0.153433683 | 1 | -0.447451588 | 2 |
| Cmma106 | 4 | -0.153436796 | 1 | -0.447450881 | 2 |
| Cmma114 | 4 | -0.153436796 | 1 | -0.447450881 | 2 |
| Cmma127 | 4 | -0.153433683 | 1 | -0.447451588 | 2 |

The small references also leave raw negative directions in the atomic-only
release blocks. Merely freezing the open cell does not cure those omissions.
This is a local partition diagnostic, not a statement that every model based
on that parent phase fails. In particular, the strongest complete optical T
reference leaves only three open-cell directions, which are raw positive at
both centers. Its zero-dimensional atomic release is explicitly reported as
empty, not as a material stability validation.

For that raw-positive joint block only, direct Schur algebra retains **two**
negative directions in the33-dimensional retained matrix. Thus an admissible
release partition does not turn the ordinary highest image into a certified
index-one TS. The centers are nonstationary; neither index, released branch
nor TS is certified. Weak-reference omission and a nonstationary/multiple-
negative center are separate issues. Do not conflate either with a failed
DFT run, or tighten the ordinary.10NEB threshold to redefine acceptance.

## Fairness and scientific boundary

Rank4T/Cmma comparisons and the full33-direction strong control remain
distinct. This is not a same-dimension prediction-accuracy benchmark. The
33-direction reference has the same complete measured data available here;
no preparation matrix or source preparation is declared free. Different
references partition that common matrix, not independent DFT replicates.
The two pairs of similar/equal Cmma scalar results are retained individually,
not converted into statistical uncertainty bars.

The geometric coverage of the two raw negative directions is exported only
as a projection diagnostic, not an energy contribution or prediction error.
Stable-R Schur algebra is never evaluated on a raw nonpositive release block;
no pseudoinverse, eigenvalue clipping or silent dimension promotion occurs.
Conversely, failure of a chosen Schur partition is not failure of an available
nonsingular direct full-internal harmonic reference. Complete harmonic
stationarity at common H/g/B remains reference-partition invariant.

There is no demonstrated adaptive-model accuracy advantage in this result.
The next physical decisions still require the second amplitude, finite-step/
energy-gradient checks and actual finite conditional branches. The registered
two-center/two-face budget, all electronic inputs and sealed physical+.005
condition remain unchanged; no final axes, domains or forecasts are selected.

## Code and provenance

`scripts/analyze_hfo2_G3_reference_complements.py` validates the E155 native
prefix reports and E150 manifest, then reuses the same pinned reference bank
as E154. `load_registered_reference_bank` is a behavior-preserving extraction
from the previous E154 verifier, committed separately as a5bd18e. All14E154
scalar and source-input fields replay identically; only that helper's source
SHA changes. The22existing adapter checks pass without modification.

Original source data remain private, read-only inputs at the recorded author
source hashes. No raw or transformed author arrays are exported. The scalar
report records code/source identities and the two exact native-report hashes:
T1cc4269a8f1437acddb94c9a7c264f52bae888669db2bcc36521982d280666ea;
preserving9bfeb3c8a7eca7533246ffd61e6e43df17262fb8ef274dfa87965c1301a10efa.
The offline helper does not claim another raw-SCF audit.

Example, replacing PRIVATE_SOURCE and FRESH_OUTPUT with local paths:

```text
python -m scripts.analyze_hfo2_G3_reference_complements
  --source-root PRIVATE_SOURCE
  --T-report benchmarks/hfo2_channels/20261008/G3_first_amplitude_E155_20261014/T_native.json
  --T-sha256 1cc4269a8f1437acddb94c9a7c264f52bae888669db2bcc36521982d280666ea
  --preserving-report benchmarks/hfo2_channels/20261008/G3_first_amplitude_E155_20261014/preserving_native.json
  --preserving-sha256 9bfeb3c8a7eca7533246ffd61e6e43df17262fb8ef274dfa87965c1301a10efa
  --output FRESH_OUTPUT
```

The first56-check test attempt had55passes/1failure: a test fixture returned
a retained-span view of its measured-basis array, so corrupting the test span
also changed the supposed fixed measurement. Return a copy in that fixture;
production arithmetic/source data were not altered. Original report retained.
Final local and independent scope, actual private-input CLI replay and source
identities are recorded in validation.json, not claimed as full regression.
