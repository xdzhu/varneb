# GaN 45.7-GPa VASP center stationarity decision loop

## Question and immutable contract

Test, rather than assume, whether the shared fixed image of the converged
B4→candidate→B1 VCNEB chain is a certifiable variable-cell stationary
index-one saddle. The common static is job 27792422. Preserve VASP 6.3.2,
PBE/Ga_d+N PAW, ENCUT 600 eV, Γ-centered 8×8×6, EDIFF 1e-7, ISYM=-1,
SYMPREC=1e-4, and external pressure 45.7 GPa. Do not tune those inputs to
make a certificate pass. Licensed POTCAR stays on hf.

## Baseline evidence and measurement model

The two ordinary VCNEB halves are converged below 0.10 eV/Å and meet at an
identical highest sampled image. At that center, the VASP stress/force-based
translation-free gradient norm is 0.002449 eV/Å and the 15D local joint
Hessian has one negative eigenvalue (-4.033 eV/Å²). The full energy/force
check nevertheless finds approximately 0.0205–0.0238 eV/Å persistent
normal-strain secant-minus-integrated-gradient discrepancies at 0.02,
0.01 and 0.005 Å. In the joint-coordinate convention,

`dH/dq_i = V(σ_ii + P_ext)/cell_scale`, where `q_i` is in Å and
`cell_scale = 3.3982714330050063 Å`. The observed residual is about
0.30–0.35 GPa (3.0–3.5 kbar) equivalent; it exceeds the previously accepted
2-kbar pressure gate. Free energy TOTEN and extrapolated energy are equal
at the printed 1e-8-eV precision for the checked ±0.005-Å statics, so a
smearing-energy choice does not explain the mismatch. Changed plane-wave
membership at fixed cutoff is a plausible, unproven contributor. Thus the
current center is not yet a strict energy-and-stress-consistent certificate.

## Bounded experiment 1 — nearby energy-informed trial

Use the 0.005-Å measured energy slopes in the three normal-strain
coordinates and the center force/stress gradients elsewhere. Project this
hybrid gradient onto the 15 translation-free joint modes and take one
Hessian-Newton correction within a 0.004-Å total trust radius. This is a
diagnostic candidate, not an automatic replacement of the archived path.
Calculate exactly nine new same-contract statics: corrected center; six
normal-strain pairs at ±0.005 Å; and a pair at ±0.01 Å along the old unstable
mode. Reuse no altered electronic parameters. Run each static with 32 MPI
tasks on hfacnormal01 via a 9-element Slurm array capped at three concurrent
elements. Do not submit twice if the array is alive or completed.

Predeclared gates (all must pass before any stronger label): complete SCF and
energy/forces/stress at every case, raw input hashes fixed; corrected-center
max atomic force ≤0.01 eV/Å; max normal/hydrostatic VASP stress residual
≤2 kbar; energy-derived normal-pressure residual ≤2 kbar; energy-secant
versus integrated VASP gradient discrepancy ≤2 kbar; and negative unstable
mode energy curvature. These are numerical-acceptance gates, not an assertion
of exact mathematical stationarity. If the gates pass, only then consider a
fresh full 15D joint Hessian plus signed basin checks at that corrected
center. If they fail, record the failure and do not reframe a proxy as a
certificate or iterate endlessly. The separate 0.10-eV/Å NEB path-force gate
does not replace a stationary-point test.

## Reproducibility and decision log

Preparation source: `scripts/prepare_gan_600eV_ts_stationarity_trial.py`.
Execution template: `cluster/hf_gan_600eV_ts_stationarity_trial.slurm`.
Raw-audit and decision results will be appended to `autoresearch.jsonl` and
summarized here, with original OUTCAR hashes and units. A normal Slurm wait
is not a failure. No CI or 235/cu17 jobs. Push only validated, selective
stage artifacts without sweeping unrelated worktree changes.

## Result and decision (2026-10-04)

Slurm array 27851541 completed all nine statics (0 exit codes). The raw
auditor verified the fixed INCAR/KPOINTS/POTCAR hashes, static geometries,
electronic convergence, output footers, and finite energy/force/stress for
every case. `manifest.json` records the proposed 0.00305864-Å joint-coordinate
Newton displacement; `stationarity_audit.json` records all nine OUTCAR hashes.

At the corrected center, max atomic force is 0.0000666 eV/Å. Independent
±0.005-Å enthalpy secants of the three normal-strain axes give a maximum
energy-derived pressure residual of 0.143 kbar. The old unstable-mode
direction has negative energy curvature, -4.310 eV/Å², with a 0.000247-eV/Å
two-sided slope. These results make the corrected center a promising
*energy-stationary, negative-curvature candidate* in the tested directions.
They do not establish the complete Hessian index at the shifted center or
the two signed basin connections.

The VASP stress residual is 3.547 kbar, and the maximum energy-secant versus
integrated-stress mismatch is 3.487 kbar. Both exceed the 2-kbar gates fixed
before submission. **Decision: the shared two-segment center is not a
certified variable-cell stationary point/first-order saddle under the
unchanged production contract.** Do not silently substitute the nearby
trial center into the already converged VCNEB chain or advertise it as a TS.
Per the bounded protocol, no fresh Hessian or further adaptive correction
was launched after this failed gate. The 0.10-eV/Å NEB convergence and
reported measured barriers remain path-level results, not a TS certificate.

The nearly unchanged 3--3.5-kbar mismatch across the old 0.02/0.01/0.005-Å
steps and this shifted-center 0.005-Å step argues against simple finite-step
truncation. VASP documents basis-set Pulay stress and basis-size sensitivity
in variable-cell calculations; changed plane-wave membership is a plausible
contributor here, not a proven unique cause. The two-sided stress/energy
check exposes a limitation of this numerical contract. Any later strict
certificate would need a separately registered strategy to resolve this
inconsistency while preserving a coherent energy/gradient definition; it
must not be manufactured by retroactively loosening the 2-kbar gate or
arbitrarily changing the 600-eV cutoff.

Relevant official documentation: [VASP Pulay stress](https://vasp.at/wiki/Pulay_stress)
and [VASP volume relaxation](https://vasp.at/wiki/index.php/Volume_relaxation).
