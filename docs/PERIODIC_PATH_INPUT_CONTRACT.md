# Ordered periodic paths: an input contract

VCNEB uses unwrapped atomic coordinates in a declared reference-cell metric.
Independently wrapping each image into its unit cell preserves its periodic
DFT geometry but can insert false lattice-vector jumps into the path. These
jumps change tangents and spring forces; correct energy/force/stress inputs
alone cannot prevent this geometrical failure.

`--mic` chooses minimum-image *endpoint interpolation*. It does not silently
modify an `--initial-chain` or `--resume-snapshot`. For a supplied chain whose
adjacent steps are explicitly intended to be short, use
`--require-continuous-periodic-lift`. This guard rejects inconsistent lifts
before loading any calculator; it never repairs the path while computing.

Prepare a short-step lift explicitly, before submission:

```python
from ase.io import read, write
from vcneb import minimum_image_path_lift, validate_periodic_path_lift

original = read("ordered_seed.traj", index=":")
images, audit = minimum_image_path_lift(original)
validate_periodic_path_lift(images)
write("continuous_seed.traj", images)
# Save audit with the immutable run inputs and retain original separately.
```

The helper preserves ordered species, cells, PBC and each physical periodic
geometry. It only adds recorded per-atom integer lattice translations and
detaches calculators; physical equivalence is not permission to use unaudited
old results. Cache reuse separately checks the exact ordered periodic geometry,
physical input hashes and complete raw SCF. It does not infer atom permutations,
rotate cells, symmetrize structures or change a Hamiltonian.

The convention is component-wise adjacent fractional minimum image, not
Cartesian nearest-image search in arbitrarily skewed cells. Half-cell steps
are ambiguous and rejected. A resolved winding is retained by following
successive images; a sparse path with unknown winding cannot be reconstructed
from endpoint equivalence. Such a path needs explicit branch information or
denser supplied geometries, not automatic nearest-atom remapping.

Always compare the segment lengths actually used by the optimizer with those
used for plotting. A plotting routine that independently unwraps coordinates
can hide a broken input band. The HfO₂ case preparation scripts now make this
comparison, and their Slurm pilot uses the guard. Generic historical defaults
are otherwise unchanged, so intentional long-step lifts are not silently
reinterpreted.

For a stopped run with a defective lift, preserve all outputs, explicitly
select a complete snapshot, register the integer shifts, reattach only exact
audited SCFs and reset the optimizer state. Do not call this an acceleration
algorithm or claim its old residual as physical convergence evidence.
