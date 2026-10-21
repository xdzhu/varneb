# Final HfO2 strain–pattern cross comparison

This is the last author-authorized two-condition experiment, specified in
[protocol.md](protocol.md). Existing minus-sign joint-release candidates
supply the two diagonal conditions; only the two missing crosses are new.
`prepared/config.json` binds their native sources and unchanged physical
contract. The prescribed physical four-vector is transported, not simply
copied as Cartesian displacements between different reference cells.

The two independently checked start norms are 0.671081 and 0.859427 A,
inside the existing 1.0-A numerical window. The four observable shifts are
opposites, predominantly Qx +/-0.246490 A; the other three are also retained
exactly, not rounded to zero. The common original-T metric independently
checks the transferred vector. The source sign was fixed before calculation.

Only `worker.py`/`policy.py` and a copied native-cache reader are needed on
the immutable E268 runtime. No runtime source or original result is edited.
Prepare and test locally before packaging; deploy to the unique E277 root,
then run native zero-SCF preflight and held-resource checks. Two 1-node,
32-MPI, 2-hour hfacnormal01 allocations are permitted, 40 fresh SCFs each.

The two original jobs completed on October 21 at 13:12:30 and 13:01:30 CST.
Only 10+4 new SCFs were used, totaling 15.52 allocation core-hours. Both
crosses met the unchanged 0.03 eV/A local target. The terminal archive
contains 114 hashed files, including every new native E/F/stress result;
all fourteen were independently replayed locally. Two diagonal J records
were also replayed, not recalculated. This finite calculation is closed.

| Strain | Held vector | Energy above same-strain well (meV/HfO2) | Released gradient (eV/A) | New SCFs |
| --- | --- | ---: | ---: | ---: |
| 0 | q0 | 83.211101 | 0.011617 | 0 |
| 0 | q1 | 114.553122 | 0.015790 | 10 |
| +1% | q0 | 92.398077 | 0.023149 | 4 |
| +1% | q1 | 107.767116 | 0.022660 | 0 |

The fixed-pattern strain responses are +9.186976 (q0) and -6.786007 (q1)
meV/HfO2. The same-strain pattern contrasts are 31.342021 and 15.369039;
their interaction is **-15.972982 meV/HfO2**, independent of initial-well
shifts. The +1% q0 cross is only 0.699140 above the existing free-chain
sampled peak; these local candidates are not new path barriers. Changing
the four-vector changes the permitted transverse and open-axis response;
this is not a unique single-mode or atomic/elastic energy partition.

## Offline reproduction

For a standard-library-only export from the source checkout:

```sh
python -m examples.hfo2_mechanism_results --output /fresh/mechanism_tables
```

For complete native replay and the main figure (existing project scientific
dependencies required), use fresh output directories:

```sh
python benchmarks/hfo2_channels/20261008/strain_pattern_cross_E277_20261021/analyze.py --output /fresh/analysis
python benchmarks/hfo2_channels/20261008/strain_pattern_cross_E277_20261021/plot.py --source /fresh/analysis/analysis.json --output /fresh/figures
```

Neither command starts DFT, modifies a scheduler, or changes the runtime.
`analysis/analysis.json` and `figures/` contain the completed results.
The complete protocol and original immutable package remain the source of
execution details, not instructions to resubmit the finished jobs.

## Historical preparation

The following commands only prepare or test and require fresh outputs:

```sh
python benchmarks/hfo2_channels/20261008/strain_pattern_cross_E277_20261021/prepare.py --output /fresh/prepared
python -m pytest -q tests/test_hfo2_strain_pattern_cross.py
python benchmarks/hfo2_channels/20261008/strain_pattern_cross_E277_20261021/package.py --output /fresh/package
```

An initial local preparation attempt exposed the projection API's requirement
for two images; single-geometry projection now uses two identical copies.
No DFT was involved. A test originally demanded bit equality for fixed rows;
the chart's round-trip differences were 1.5e-36 A and are tested with the
existing 1e-12-A numerical tolerance, without changing geometry or physics.

The two terminal outcomes, not a desired sign or journal requirement, end
this calculation protocol. Missing candidates remain missing; no automatic
expansion is allowed. Old paths and independent-prediction labels remain.
