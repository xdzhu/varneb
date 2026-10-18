# E274: four frozen midpoint atomic/cell crossovers

**Completed28950839, C0 at2026-10-19 03:13:44 CST onnode1.** Exactly four
new frozen SCFs used6m33s, 3.466704native-call proxy/3.493333allocated
core-hours. The original AA/BB caches were reused. Sourcece3e056d was
pushed before deployment; the16-file immutable bundle pin is
`d9cc4c2b81ecf7dabb542c943913a75783e5e0a1d472c0a47c591588bbb4f6ca`.
The [finite protocol](protocol.md) ended at these four calls without retry,
relaxation, new chain or larger matrix. The complete36-file native archive
SHA256 is`cb9aebee004d2ed744acf0e941b85b114040b9913fad012dde0a490b6d4e1446`.

## Result

| Contrast, meV/HfO2 | Preserving | PO-to-T |
| --- | ---: | ---: |
| BB minus AA | 10.820844 | -0.131935 |
| Internal change at H_A | 11.398469 | -0.056351 |
| Internal change at H_B | 10.362204 | -0.081225 |
| Cell change at s_A | 0.458639 | -0.050711 |
| Cell change at s_B | -0.577625 | -0.075584 |
| Finite interaction | -1.036265 | -0.024874 |

The preserving initialization difference persists at both fixed cells.
Changing the open row alone has a much smaller, pattern-dependent effect.
This rejects an open-cell-only explanation of that finite gap, not every
elastic role along the path. Internal and cell interventions do not commute
as additive energy components; their interaction is not a fitted mode
coupling coefficient. T's conditional effects remain small under this test.
The frozen contrast still does not distinguish a separate continuous branch
from residual atomic relaxation or supply a newly relaxed barrier/TS.

The [analysis](analysis/analysis.json) independently replays every new
native E/F/stress, exact native input proof and four diagonal energies.
Pseudo/orbital bytes were checked on HF, not redistributed or rehashed
locally. Source tables retain all eight corners, marking which four are new.
Older path, endpoint, curvature and release costs remain additional;
no old forecast label or score is replaced.

```text
python benchmarks/hfo2_channels/20261008/midpoint_factorial_E274_20261019/analyze_terminal.py --output FRESH_OUTPUT_DIRECTORY
```

## Admission and original preparation

Prepared from the **completed** E268/E269 sampled peaks, not the earlier
unconverged snapshots. Exactly four new static SCFs were admitted.

The preserving starting-band contrast is10.820844meV/HfO2; T's is
-0.131935. The open-cell row and internal pattern both change. Exchange
them in the same ordered fractional chart, keeping the two existing AA/BB
energies. AB/BA are frozen geometries, not relaxed candidates. This
distinguishes conditional internal and affine-cell effects without a
unique phonon/elastic energy partition or a new independent forecast.

All four native INPUT/KPT bytes remain identical to the original contract.
No pseudo/orbital files are distributed. The geometry preparer checks the
same substrate, ordered fractional registration, minimum separation and
native STRU roundtrip. Actual HF/ASE preflight, held resource validation
and accepted-handle release were separate required gates; their receipts
are retained here. The first native SCF took95.303646s with genuine32MPI,
unchanged six input files and complete E/F/stress.
The native writer's periodic wrapping is explicitly distinguished from
the continuous analysis lift. Formal calls must match the STRU bytes
registered by the native preflight, not an earlier serialization's hash.
This was caught before any submission or DFT call, without changing physics.

One32CPU hfacnormal01 allocation, one hour, serial mpirun-np32 and one
thread/rank; no CI, new NEB, relaxations, retry, model fit or extended grid.
Run code is frozen in a new namespace; prior runtimes are read-only.
Original labels, unsuccessful releases and preparation costs remain.

`control.py` implements native preflight, duplicate-resistant held dispatch,
the four-call compute-only entry and C0-only portable collection. The
existing package helper copies only this control, run script, protocol and
prepared decks into a byte-pinned bundle. No running source is replaced.
