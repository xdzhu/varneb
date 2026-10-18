# Strain response of the lower connected preserving candidate

**+0.5% complete; +1% still running.** Original28951226 completed0:0 at
04:06:49 with7fresh interior SCFs, zero optimizer updates and ordinary
residual0.0923153015eV/A. The native sampled peak is87.7185447330meV/fu,
10.8123321834 below the registered98.5308769164. T escape is90.4181783358,
so the new sampled separation is2.6996336028, not a resolved physical ordering.
The transverse structural proxy remains0.1210822924A. An initial force pass
does not prove initialization-independent relaxation or branch stability.

All7fresh calls and the selected9native E/F/stress records have been
independently replayed locally. The128-file archive is285005bytes, SHA256
`0eec746dc93431d001a57f72ebea3f38f941fbd62b514fa005d9684d475de02b`.
New cost6.128975native proxy/6.177778allocation core-hours; source path,
release and earlier strain calculations remain additional. Full assets
were checked onHF, not redistributed or rehashed locally.

The [four-panel mechanism comparison](figures/strain_0050/README.md)
retains all four same-strain candidates, including the prior upper-seed
control. `completed_analysis_0050` explicitly marks+1% pending and leaves
old labels/scores untouched. It is not a new independent prediction.
The original+1% handle28951227 continues with no retry or runtime change.

The admission/startup records below remain historical evidence.

Actual jobs28951226(+0.5%) and
28951227(+1%) started2026-10-19 at03:55:14CST onnode1/node13. Each is
the verified one-node32CPU/18h allocation, released only after native
zero-SCF preflight and held identity/resource checks. Runtime source
4e1fdf26 was pushed before deployment; bundle pin
`00ea6f9f82e50aeea60cfd560ecdcc76cec68b3308a0683e8bdf0fe922dfe8e3`.
No prior runtime package was overwritten.

At03:58:55 each had two completed fresh SCFs and no failure record. The
first calls took98.111/98.864seconds; fixed six assets, real32MPI and
native E/F/stress replay passed. This is a startup prefix, not a completed
band or convergence claim. See accepted/held/released, preflight,
package and startup receipts. The read-only inspector lives separately
from the immutable runtime bundle onHF.

The finite [protocol](protocol.md) transfers the actual E273 step9 lower
zero-strain candidate to the already studied +0.5% and +1% boundaries,
reusing their original fixed endpoint caches. The two complete free bands
ask whether transverse atomic accommodation changes the strain-dependent
competition with T escape. This is not another forecast competition.

Two one-node32CPU allocations at most, ordinary0.10/noCI,60updates/18h
per band,427fresh interior calls maximum per band. All original electronic
settings and data labels remain unchanged. No automatic extra study follows.
The scheduler's walltime is an upper execution limit, not a convergence ETA.

```text
python benchmarks/hfo2_channels/20261008/preserving_strain_E275_20261019/prepare.py --output FRESH_DIRECTORY
python -m pytest -q tests/test_hfo2_preserving_strain.py
```
