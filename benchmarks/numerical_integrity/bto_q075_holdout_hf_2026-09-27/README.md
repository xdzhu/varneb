# BTO transverse-soft independent holdout `(Q_z,Q_x)=(0.75,0.15)`

This coordinate was specified in the integrated submission plan **before**
the holdout energy was computed. It is independent of the four fitted
conditional candidates and is intended to test interpolation and branch
selection, not to fill in a contour by assumption.

The three-start preflight fixes the five-atom `1×1×1` cubic Γ-mode source,
ABACUS PBE/100 Ry/Ba–Ti–O 10 au DZP, `4×4×4` **electronic** k mesh, and
the `Q_y=0,+,-` seeds. Slurm `27787471` on hf `hfacnormal01` completed the
three static seeds with 32 MPI; the independent raw INPUT/KPT/STRU/SCF/
force/stress audit passed. The signed seeds are degenerate to numerical
precision and each lies `12.8465 meV/BTO` below the frozen seed. This is
only a static branch-start observation.

After a calculator-free resume preflight, Slurm `27787487` began the
three-start conditional relaxation with the same calculator contract,
orthogonal gradient target `0.003 eV/(sqrt(amu) Å)`, and `1 kbar` center
stress target. Its result is **pending**. A completed Slurm job alone will
not certify stationarity, curvature, branch identity, or interpolated PES.
The next gate is a full raw-output audit and replay of all three branch
outcomes, followed by a comparison of this held-out energy with an
interpolation trained without this point.

The JSON files in this folder are byte-for-byte copies of the signed source
under `/public/home/iai806/abacus/agent-runs/20260927-varneb-bto-soft-gridpilot-qx030`;
the large raw calculator tree remains on hf. The exact frozen production
source differences are archived in `../bto_hf_source_snapshot_2026-09-27/`.
No HfO₂, GaN, or additional
BTO protocol is inferred from this holdout.
