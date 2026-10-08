# Exploratory reuse of the existing G0 joint probes

This analysis reuses the eight static calculations of job28240257. It does
not add an SCF, alter a production input or substitute for G2 clamped training.
Its center is the historical guided-chain image02, which is nonstationary.
The first direction is a normalized atomic chain secant, **not a phonon**;
the second is scaled symmetric xx+yy strain. Both use the original free-cell
chart and cell scale5.127945044745648 Angstrom, at zero pressure.

The previously archived diagonal-curvature check is retained unchanged. The
new analysis assembles both complete42-coordinate gradient differences and
their mixed projection, using the public directional joint-curvature API.

```console
python -m scripts.analyze_hfo2_joint_probe_reuse --output /new/path/local.json
python -m scripts.analyze_hfo2_joint_probe_reuse --raw-root /historical/probes --output /new/path/raw-audit.json
```

The raw mode rechecks the source/center log and geometry, byte-identical
INPUT/KPT/pseudopotentials/orbitals, every signed chart geometry, genuine32MPI
SCF convergence, and full energy/forces/stress. Its rebuilt gradients and log
hashes must agree with the original eight archived point audits. It writes
only a new report, never the source calculations or old audits.

Manifest JSON content must agree; both actual raw and archive manifest hashes
are recorded because Git can normalize historical JSON line endings. This
does not relax the mandatory raw-byte hashes of any DFT input or source log.

At0.01 Angstrom the raw two-direction projection is, in eV/Angstrom2:

```text
 4.671228720   -3.318399815
-3.318981707   19.070358556
```

At0.02 Angstrom the mixed entries are-3.317902953/-3.319421332. The raw
reciprocity defects are4.07642e-5/1.06942e-4. Two-step full-action operator
spread0.117872 eV/Angstrom2 is retained as a diagnostic, not a total error bar.

Positive projected eigenvalues(3.94316,19.79843) do not prove full stability.
Transverse column norms5.86712/7.09425 eV/Angstrom2 demonstrate substantial
coupling outside this slice; their self-curvatures remain unmeasured. The
center gradient norm0.118695 eV/Angstrom also precludes stationarity claims.
No omitted block is eliminated, no conditional surface is asserted, and this
post hoc reuse is not an independent prospective barrier prediction.

The [local archive analysis](local_clean_analysis.json) and
[actual HF raw-source analysis](hf_raw_analysis.json) are retained separately.
The [delivery receipt](validation_delivery.json) records903passed/2skipped
in88.61s, fresh center+8SCF checks, all eight original log hashes, and identical
local/HF full Hessian actions. The same-archive analysis/API source hashes agree.
The later converged-band dataset has its own regression/delivery record.
