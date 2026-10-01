# GaN original B4 endpoint, isolated strict native VASP continuation

This is the archived result of Slurm job `27812112` on `hf` (completed 0:0,
2026-10-01). Source workspace:
`/public/home/iai806/abacus/agent-runs/20261001-varneb-gan-b4-endpoint-600eV/run2`.
The input POSCAR is the unchanged CONTCAR of the prior one-frame native test
(job `27812043`), itself the original production-chain B4 image 0.

The physical electronic contract remains VASP 6.3.2, PBE, Ga_d+N PAW,
600 eV, Γ-centered 8×8×6, `EDIFF=1e-7`, `ISYM=-1`, `SYMPREC=1e-4`, and
`PSTRESS=457` kbar. Relative to the first native test, only the ionic stop
changed (`EDIFFG=-0.02` to `-0.005` eV/Å) and `NSW` became 30. The POTCAR
is intentionally not redistributed; its SHA-256 is recorded in the manifest.

The raw OUTCAR has five evaluated ionic frames and clean SCF/run-end markers.
The final maximum atomic force is 0.000151 eV/Å and the maximum raw-stress
residual from 45.7 GPa is 0.05888 kbar. The final B4 has at most 0.00225 Å
relative atomic displacement and 4.11×10⁻⁵ relative volume change against
original production image 0. Its native enthalpy shift is -0.03598 meV/GaN.
See `benchmarks/numerical_integrity/gan_45p7_vasp_b4_endpoint_strict_run2_20261001.json`
for input/output hashes and the executable audit.

This is an endpoint diagnostic, **not** a replacement for the archived VCNEB
image 0, a recalculated barrier, or a strict transition-state certificate.
