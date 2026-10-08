# Same-subspace atomic--cell curvature preparation (G3 material probes pending)

The new `vcneb.ActiveJointCurvatureCoordinates.for_clamped_plane` shares the
exact open deformation subspace already used by endpoint BFGS and VCNEB. Its
one/three cell directions are not replaced by six symmetric strain directions.
General oblique/rotated substrate vectors stay fixed in every generated probe.
All atoms remain mobile; translations are projected explicitly. Twelve-atom
tilt-released starters have 39 joint coordinates and 36 translation-free ones.

The legacy `vcneb/joint_curvature.py` is byte-identical, SHA256
`6c83d2ecbc216e56cede00b58945645807688a993ba0bce4a40b85a3dd8730c1`.
An initial attempt to extend it correctly failed the published GaN provenance
gate. The independent subclass avoids relabelling historical source/data;
neither the old manifests nor numeric results were changed. Regression also
compares default new/old structures and physical gradients exactly.

Reproduce without DFT or submission, in a fresh report path:

```bash
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
python -m scripts.verify_clamped_joint_coordinates --report fresh_verification.json
```

Verification includes all 10 checksum-validated, unrelaxed HfO2 G2 starters,
five probe geometries each, and 8 oblique/rotated/currently-deformed Cu/EMT
energy--force/stress work checks, both normal-only and tilt-released boundaries.
The 232 EMT evaluations are cheap implementation checks, **not new HfO2 DFT**.
The prespecified 3e-5 eV/A check is not a barrier-error or material-convergence
tolerance. Nonzero pressure in the Cu check does not alter the HfO2 P=0 contract.

`local_verification_final.json` records the Windows working-tree check.
Git normalizes line endings in four of its source files: use
`local_clean_verification.json` and `hf_verification.json` for the identical
exported source archive. All five module hashes match between those receipts;
the legacy hash above is identical in the working tree and both archives.
Clean regression passed **864 tests, 2 skipped, in 142.06 s**. Local ASE 3.28.0
and HF ASE 3.23.1b1 both checked 10 starters and 232 EMT evaluations. Maximum
work errors are 1.100e-9 and 1.840e-9 eV/A, respectively; HF substrate drift
is below 8.1e-20 A. See `validation_delivery.json` for the exact tree/archive.

This evidence does not yet supply a HfO2
joint Hessian, stable elimination at a bottleneck, TS certificate, favorable
strain window or independent barrier prediction. Those retain the G1/G2/G3
material gates in the research plan; no additional production job is submitted.
