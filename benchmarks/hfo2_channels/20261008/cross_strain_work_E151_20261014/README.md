# E151: observed cross-strain work, not a stationary barrier derivative

This zero-DFT analysis uses all eight admitted0/+1% training bands,72frame
records/66distinct native log-STRU identities and36equal-image-index pairs.
The exact raw initial cache is common to all four candidates at each strain.
All structures use one ordered unwrapped PO+ reference chart, with the original
fixed-plane/open-tilt space,5.12968067458423Angstrom cell/control scales and
raw physical E/F/stress. No source ASR, mode/reference selection or extra SCF.
The admission audit and pairing reader each reparse72records; these are not
144new calculations. Undistributed full pseudo/orbital bytes are not locally
reaudited by this reader; their historical native audit remains the source.

## Defined numerical connection

For a matched pair, delta_y=y(+1%)-y(0) and mean_g=(g0+g1)/2 in the same
dual chart. The endpoint trapezoid is mean_g dot delta_y. Split it into
internal atomic, explicit rigid translation, permitted open row and prescribed
substrate terms. Subtract the same common initial-pair work when comparing a
sampled peak. Independent Cartesian force/raw-stress contractions agree with
the dual-coordinate contractions within6.662e-16eV/cell.

No intermediate structure is generated: the connection is a declared
straight fractional-atom/open-row/strain connection, not a reoptimized MEP
in strain. No+0.005geometry/label is opened. All paired native peak indices
match between the conditions(M3; T/each flip4), so their actual energy
differences reproduce the observed sampled barrier shifts. Equal image index
is not an automatically canonical continuous MEP correspondence.

## Actual matched-peak contrasts

All columns are meV/HfO2 relative to the common initial-pair change. Atomic,
open-row and substrate columns are connection-dependent signed work, not
state-function mode energies. Translation terms are retained in the data
(maximum peak magnitude2.411e-10), though negligible at this display precision.

| Candidate | Observed sampled peak change | Atomic work | Open-row work | Substrate work | Energy minus full work |
|---|---:|---:|---:|---:|---:|
| PO-to-M | -1.230403 | 2.720687 | -0.174328 | -3.835247 | 0.058485 |
| PO-to-T | 12.347027 | -0.228947 | 0.220068 | 12.275095 | 0.080811 |
| Preserving flip | 28.596248 | 36.895288 | -0.078916 | -12.858888 | 4.638764 |
| Reversing flip | -44.024331 | -32.332083 | -0.270001 | 0.072814 | -11.495061 |

The prescribed-coordinate term dominates the T candidate's declared
connection. For preserving, its sign is opposite to the observed change,
with a substantial atomic term and nonzero quadrature residual. Thus an
isolated mechanical partial derivative is not the observed barrier response.
Neither a small closure residual nor this work split proves a causal mode
coupling, stable branch, continuous TS, DFT/barrier error interval or H2.
The two already running E150 local Hessians/conditional tests remain needed.

## Strong-reference comparison guard

The existing full-release harmonic solver must agree with a direct fixed-
control solve under complete retained/released changes of basis. Six checks
cover retained ranks1/3/4 and rotated coordinates, nonzero gradients and an
index-one internal block; the prescribed control is excluded from that index.
These are analytic consistency checks, not material predictions or a new
theorem. At this complete harmonic limit a reference rename cannot improve
accuracy. Any claimed advantage must instead survive the registered strong
fixed-reference comparison at measured resolution/cost, including validity,
nonlinearity and branch tests. No strong baseline is dropped or weakened.

## Reproduction and actual verification

```text
python -m scripts.analyze_hfo2_cross_strain_work --output NEW_REPORT.json
python -m pytest -q tests/test_hfo2_cross_strain_work.py tests/test_stationary_reference_invariance.py tests/test_stationary_branch.py tests/test_strain_work.py
```

Actual local related tests:72pass/0failure/error/skip. Sources test native
contrasts, independent EMT physical-work identities, raw nonzero net forces,
long lifts, missing training sources and rejection of an initial integer
remapping. The first exact dictionary equality was too strict for binary64
trajectory serialization differences4.441e-16Angstrom, despite identical raw
initial caches. Actual source geometry/energy were not changed: ordered-lift
agreement and work identities now allow only explicit1e-12algebraic roundoff,
not a physical uncertainty bound; nontrivial image remapping is still refused.
Final native report is regenerated with the final admission guard. The
preliminary export is retained outside Git in E:/TEMP/varneb-E151-preliminary-analysis-20261014.json;
all scientific fields are exactly equal, and its source is reproduced by
removing only the subsequently added two-line unique-eight-source guard.
Independent execution:72pass/0failure/error/skip,13.354JUnit seconds on
tree646b7cfdc2e80448ddaa98bd6855d429b6424079; all6449tracked files were
present before testing. Archive SHA256 is
a0b88df333ea1377c475640978dbc8aea6ce22a0a094452ce2abdf48ac4c9017.
The actual independent CLI replay reproduces every scientific field exactly.
Only script_SHA metadata differs: verified LF/CRLF source export, not a
physical-data normalization. Both original exports and JUnit reports are
retained here. No whole E151 regression, new material SCF or paper edit is
claimed. See validation.json for the executed-byte hashes and scope.

At19:20:33CST, native28855614/28855615 both remain RUNNING, each23/148
completed and24started,32CPUs/node10,node27; stderr0. Completed receipt times
average102.655/104.513s. This status reads scheduler/evaluation receipts, not
a new raw-SCF reparse. The previously native-audited first evaluations remain
the physical startup evidence. Rough completion23:00–01:00CST is conditional
on comparable later cost. No scheduler writes, extra job, live-source changes
or reserved physical+0.005geometry. Full research goal remains incomplete.
