# GaN 600-eV normal-strain derivative step-size check

Slurm array `27812251_0`–`27812251_11` completed 0:0 on `hf`, 2026-10-01.
The isolated source was
`/public/home/iai806/abacus/agent-runs/20261001-varneb-gan600-normal-strain-steps/work`.
These are 12 static VASP evaluations: axes 12–14 of the established joint
atom–strain chart, at ±0.01 and ±0.005 Å about the same one-step-refined
image-15 center used in the earlier ±0.02-Å Hessian audit.

Every worker checked its input hashes before running. The calculator contract
was unchanged: VASP 6.3.2, PBE, Ga_d+N PAW, `ENCUT=600 eV`, Γ-centered 8×8×6,
`EDIFF=1e-7`, `ISYM=-1`, `SYMPREC=1e-4`, and static `IBRION=-1`, `NSW=0`.
All enthalpies use the same *analysis* pressure, 45.7 GPa. POTCAR and
WAVECAR are intentionally not redistributed. The POTCAR SHA-256 recorded in
the manifest and each case was checked on `hf`; this archive preserves the
other inputs and the raw OUTCAR for every case.

From repository root, reproduce the report with:

```bash
python -m scripts.audit_gan_600eV_normal_strain_steps \
  --output /path/to/fresh/normal_strain_audit.json
```

The committed report is
`benchmarks/numerical_integrity/gan_600eV_normal_strain_step_dependence_20261001.json`.
It compares each enthalpy secant with a Simpson-integrated stress/force
gradient. Across the three axes, the residual remains 0.0205–0.0238 eV/Å
for step sizes 0.02, 0.01, and 0.005 Å. The 0.005-Å pairs still change their
plane-wave counts at 355–361 of the 384 k points, even when the single
reported maximum plane-wave number is unchanged. The mismatch therefore does
not collapse with a smaller finite-difference step. This observation is
compatible with a variable finite basis; it neither uniquely proves a Pulay
mechanism nor certifies a stationary transition state. No VCNEB endpoint,
path, cutoff, or published barrier was replaced.
