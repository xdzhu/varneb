# E073: allowed gradients versus prescribed-coordinate reactions

The ordinary-converged zero-strain PO-to-M chain is replayed from the nine
portable native E068 INPUT/KPT/STRU/log/E/F/stress records. No SCF, optimizer,
geometry, physical parameter, model selection or holdout access occurs.
The six complete physical-file checks remain historical HF evidence, not
new checks of the proprietary bytes on this workstation.

The public saddle summary previously exposed active **raw** force norms,
whereas the per-image path tangent/perpendicular fields already used the
implemented constraint projector. Large forces excluded by the prescribed
clamped plane can therefore dominate the legacy summary without being a
stationarity failure in the actual allowed space. The production NEB force
was already projected correctly; it is unchanged.

New explicit `allowed_true_*` and `constraint_reaction_*` fields are additive.
Old fields and their numerical meanings remain for reproducibility. They
represent generalized forces under the original cell scale, not raw kbar
stress, a tolerance conversion or a Hessian measurement. The reaction is
the projector-excluded **active** force; it does not recover forces already
removed by a separate coordinate mask. This applies only to physically
enforced optimization constraints, not a post-analysis truncated mode basis.
Support reactions remain physical: they must still enter the external
substrate-work derivative through the raw stress and actual cell Jacobian,
even though they are excluded from internal stationarity diagnostics.

At sampled image3, the legacy active raw maximum is2.0650935345eV/Angstrom;
the excluded maximum is2.0631303934 and the allowed maximum is0.0900238083.
The allowed Euclidean norm is0.2317176664 and signed physical tangent is
-0.1242638810; these are different norms, not contradictory tolerances.
The ordinary max-vector residual remains0.0992140881. The allowed physical
perpendicular norm0.1955800726 differs from the raw legacy value2.6171263383.
Diagnostics reproduce the ordinary force byte-for-byte before/after replay.

These measurements prevent rejection of a clamped calculation because of
its support reactions. They do not certify an exact stationary bottleneck,
negative-mode index, continuous maximum, energy error or H1/H2 prediction.
The reported path finite-difference curvature is not a local Hessian.
No tighter ordinary NEB threshold, CI or extra phonon/Hessian budget is added.

Reproduction, from the repository root into a **new** output:

```text
python -m scripts.analyze_hfo2_clamped_stationarity --case benchmarks/hfo2_channels/20261008/clamped_M_terminal_E068_20261011 --audit-sha256 56b92e8b15195d6f619a8402c7eee83a774b2e364e012478faa3fdc5d2da431c --output NEW_REPORT.json
```

`analysis.json` preserves source/core/receipt/observation hashes and all nine
new diagnostics. Separate validation records distinguish actual material
replay from synthetic tests and current HF state. The live immutable HF
runtimes are not overwritten by this local diagnostic-interface improvement.
