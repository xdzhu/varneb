# Strain response of the lower connected preserving candidate

**Both original targets completed0:0.** Job28951226(+0.5%) ended04:06:49;
28951227(+1%) ended04:51:10 on2026-10-19 CST. All42fresh native interior
energy/force/stress calls and both selected nine-image bands have been
independently replayed locally. No restart, parameter change or endpoint
recalculation was required.

| Strain | Registered preserving | Released-source preserving | Registered T | T minus new preserving |
|---|---:|---:|---:|---:|
| 0 | 96.828183 | 83.166281 | 84.416624 | 1.250343 |
| +0.5% | 98.530877 | 87.718545 | 90.418178 | 2.699634 |
| +1% | 125.424431 | 91.698937 | 96.763650 | 5.064714 |

Energies are native sampled peaks in meV/HfO2, not certified stationary
barriers. The new0-to-+1% rise is8.532656 rather than28.596248meV/fu;
the20.063593difference shows candidate-response sensitivity to internal
accommodation. Small positive separations do not resolve global physical
ordering. One prescribed source, not a unique continuous branch or a new
independent prediction, was evaluated. Old labels and all scores remain.

The midpoint passed at step0 after7fresh calls, residual.0923153015eV/A;
the high-strain band required4updates/35calls, residual.0909313638eV/A.
The midpoint's initial pass is not proof of relaxed-branch uniqueness.
Transverse structural proxies at the peaks are.121082/.122675A, not
electronic polarization or single-phonon amplitudes. New total cost is
35.798164native-call-proxy/36.008889allocated core-hours; source path,
release and earlier strain costs remain additional and are not double-counted.
Full six assets were checked onHF, not redistributed or rehashed locally.

The complete [three-strain summary](figures/summary_three_strains/README.md)
is manuscript Figure9. The [+1% profile comparison](figures/strain_0100/README.md)
and earlier [four-candidate midpoint comparison](figures/strain_0050/README.md)
retain all sampled curves, including the upper-initialized midpoint control.
`completed_analysis_full` contains10bands/90native sample rows and both
actual target archives. `completed_analysis_0050` remains a dated partial
snapshot; its pending flag is historical, not the current job state.

The original P central samples are geometrically Pbcn at the three fixed
symmetry tolerances. The released-source zero center is Aea2; the strained
centers change Pc/Pna2_1/Aea2 with tolerance. These are post-processing
labels, not new phases or DFT symmetrization. The +1% local fixed-observable
A/J experiment remains missing despite this complete free-chain result.

Archives: +0.5%,128files/285005bytes, SHA256
`0eec746dc93431d001a57f72ebea3f38f941fbd62b514fa005d9684d475de02b`;
+1%,324files/720284bytes, SHA256
`cd6803ba0c7971d9719049c5e7dee23599fb25fb7a19cf17271e608f2186553b`.
Collection receipts preserve original Slurm rows, source/bundle hashes and
the separate collector version. Collection/replay itself creates zero SCFs.

Replay the complete data without an allocation or an ABACUS executable:

```text
python benchmarks/hfo2_channels/20261008/preserving_strain_E275_20261019/analyze_terminal.py --output FRESH_DIRECTORY
python benchmarks/hfo2_channels/20261008/preserving_strain_E275_20261019/plot_summary.py --source FRESH_DIRECTORY/analysis.json --output FRESH_FIGURE_DIRECTORY
```

The finite two-target protocol is closed. No additional local relaxation,
Hessian, strain, forecast arm or allocation is automatically scheduled.

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
