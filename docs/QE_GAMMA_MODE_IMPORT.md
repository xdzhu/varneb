# Importing a QE Gamma reference without changing the calculator contract

`vcneb.qe_modes.read_qe_gamma_modes` reads the first, explicitly zero-q,
complete `3N` block of a QE `matdyn` **flvec** file. It does not run QE, infer
a structure or finite-q units, or convert a literature calculation into an
ABACUS calculation. Supply one positive mass per atom in the original row order.

```python
from vcneb.qe_modes import read_qe_gamma_modes

reference = read_qe_gamma_modes("matdyn.modes", masses_amu)
E = reference.mass_weighted_eigenvectors  # complex (3N, 3N) columns
```

## Vector convention and numerical limits

QE flvec columns are Cartesian displacement vectors, obtained by dividing
dynamical-matrix eigenvectors by the square root of each atomic mass and
normalizing them. They need not be mutually orthogonal in a Cartesian metric.
The importer restores the mass-weighted convention by multiplying each atom's
three rows by `sqrt(mass)`, then normalizing each column. This convention is
documented in both the [QE matdyn reference](https://www.quantum-espresso.org/Doc/INPUT_MATDYN.html#flvec)
and the [official QE 6.3 source](https://raw.githubusercontent.com/QEF/q-e/qe-6.3/PHonon/PH/matdyn.f90).
This importer is **not** a reader for the distinct orthogonal fleig format.

Both printed displacements and restored vectors remain complex and read-only.
The source SHA256, signed frequencies, supplied masses, printed normalization
defect and restored Gram defect are retained. There is no implicit real-phase
choice, QR orthogonalization, acoustic sum rule, symmetry repair or force-constant
reconstruction. A common mass scale cancels during normalization; extreme,
numerically unresolved mass ratios are rejected rather than silently repaired.

Defaults accept printed normalization defect below `5e-5` and restored Gram
defect below `5e-4`; these are **file-precision checks**, not DFT/NEB tolerances.
The reader bounds input to 12 MiB and 128 atoms and rejects missing/duplicate
mode labels, malformed rows, nonfinite numbers, inconsistent THz/cm-1 labels,
and a first nonzero q. Subsequent q blocks are not imported. A syntactically
valid file is not proof of association with a structure: verify its source,
ordered atoms, masses and frame separately.

The existing exact-basis projection interface has a stricter orthogonality
contract. An accepted printed Gram defect of about `4e-6` does **not** make
the vectors an exact basis that can be passed through that interface unchanged.
An eventual real/orthonormal subspace representation needs its own explicit
phase, degeneracy and projection-error policy; do not discard imaginary parts
or silently loosen the projection contract to make a file pass.

## HfO2 reference audit

The reviewed reproduction entry point is:

```text
python -m scripts.audit_hfo2_cmma_reference --source-root AUTHOR_FILES --output NEW_REPORT.json
```

Only the five pinned files listed in the script are accepted, by Git blob
hash and byte bound; an existing output is refused. Raw author files are kept
outside Git, not bundled in VARNEB. No author script is executed. The source
of each file, proper rotations, integer basis changes, species assignment,
periodic shifts, rounding errors, symmetry tolerances and library versions
are recorded. See the [reference audit](../benchmarks/hfo2_channels/20261008/cmma_reference/README.md).

That audit currently connects two representations of the **published** Cmma
reference. Mapping it into the production ordered-atom/path gauge remains a
separate gate. Literature LDA modes can supply geometric directions, not our
PBE curvatures or barriers. This feature changes no cutoff, orbital, k mesh,
SCF setting, pressure, endpoint, production source or running job.
