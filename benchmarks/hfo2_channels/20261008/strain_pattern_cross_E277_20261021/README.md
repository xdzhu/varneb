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

At this preparation stage no new DFT result exists. The following commands
only prepare or test and require a fresh output directory:

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
