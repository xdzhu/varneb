# Directional joint-curvature implementation verification

This milestone prepares the G3 atom--cell sampling pipeline. It adds paired,
calculator-free probes and assembles measured full-chart Hessian actions from
selected directions, retaining response outside their span. It does not change
the production calculator contract or authorize additional material jobs.

Run from an isolated source archive:

```console
python -m scripts.verify_joint_curvature_stencil --report /new/path/verification.json
```

The script needs NumPy and ASE, not pytest or a DFT executable. It checks 720
geometry probes from ten registered, unrelaxed HfO2 starters and makes 180 Cu/EMT
evaluations. Its two displacement amplitudes and numerical tolerances are
implementation checks, not HfO2 production settings or DFT error estimates.

The synthetic example measures two Hessian columns while retaining coupling
to a third direction. The unmeasured third diagonal can have either sign;
therefore the projected spectrum cannot certify the full stability index.
Raw reciprocity defects remain visible before symmetrization. No unsampled
block is eliminated without its own stability evidence.

The [delivery receipt](validation_delivery.json) identifies tested tree
`a6dad26a` and archive SHA256 `45dd4d67...`: 895 tests passed, 2 skipped in
91.30 seconds from the isolated archive. The same archive's local ASE3.28.0
and HF ASE3.23.1b1 checks agree on all six module hashes. HF needed no pytest
installation and ran no DFT. Its maximum energy/gradient curvature difference
is 9.1720e-5 eV/Angstrom2, with substrate drift4.8736e-20 Angstrom.
See [local verification](local_clean_verification.json) and
[HF verification](hf_verification.json). Worktree byte hashes
can differ through Git's Windows newline normalization; compare local clean
archive and HF module hashes. The legacy joint-curvature module is unchanged.

Material Hessian measurement, transition-state certification and independent
barrier prediction are explicitly **not** established by these checks.

See [the API manual](../../../../docs/JOINT_CURVATURE_PROBES.md) for the metric,
pairing, physical-gradient and provenance requirements.
