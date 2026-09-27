# BTO off-axis conditional-branch evidence, `Q_z=0.9`

These thirty-six JSON files are byte-for-byte copies of the completed `27786180`,
`27787086`, `27787309`, and `27787447` ABACUS runs and subsequent read-only
audits on hf. The eighteen newly imported files were checked against their
remote SHA-256 hashes.
Their source directory is
`/public/home/iai806/abacus/agent-runs/20260927-varneb-bto-soft-gridpilot-qx030`.
The `run-q090_q000/` and `run-q090_q030/` subdirectories there retain the
original `INPUT`, `KPT`, `STRU`, ABACUS logs, and immutable evaluation cache.
The local copies here retain the SHA-256 links back to those raw inputs and
logs; they do not duplicate the large raw calculation tree.
The byte-identical hf production source files that differ from current
development are archived in `../bto_hf_source_snapshot_2026-09-27/`.

The two selected `-Q_y` branches satisfy the stated orthogonal-gradient and
1-kbar residual-stress gates. Independent raw-output audits cover 77 and 78
DFT evaluations. Three-start cache replay resolves both the higher-energy
frozen branch and the nearly degenerate `+Q_y`/`-Q_y` pair. The `0.05`
`sqrt(amu) Å` curvature-preflight files certify that all 32 signed probe
geometries are safe. The subsequent two array elements of `27787086` completed;
independent raw-output audits now cover all 109 and 110 cached points,
respectively. The reconstructed 16-dimensional force-difference Hessians
have no negative eigenvalues at this **one** step, with lowest values
`0.00569877` and `0.00581667 eV/(amu Å²)`. The largest energy/force diagonal
curvature discrepancies are `0.007674` and `0.010632 eV/(amu Å²)`, larger than
the lowest eigenvalues. Thus this is a useful positive screen, **not** a
numerical certificate of local-minimum stability.
The 0.10-step array `27787309` completed with `0:0` for both elements. The
independent raw-output audits now cover **141 and 142** cached DFT points.
The new 16-dimensional Hessians again have no negative eigenvalues, with
lowest values `0.00572021` and `0.00584450 eV/(amu Å²)`. Relative to 0.05,
their lowest-direction absolute overlaps are `0.99999984` and `0.99999965`.
The 0.10 energy/force diagonal curvature discrepancies fall to `0.001760`
and `0.003656 eV/(amu Å²)`; the older 0.05 discrepancies remain larger than
the lowest curvature. The old 0.05 point sets were reconstructed without new
DFT through symlink-only audited-cache views, because the append-only run
directories now also contain the 0.10 points. No original cache was changed.
These comparisons support step stability, **not** a rigorous eigenvalue-error
bound or local-minimum certificate.

The two four-probe lowest-mixed-direction geometry preflights passed without
DFT. Array `27787447` then completed at the unchanged calculator contract.
Independent INPUT/KPT/STRU/SCF/force/stress audits cover **145 and 146** total
cached points and verify each of the eight new signed probes. At steps
`0.05/0.10 sqrt(amu) Å`, respectively, the direct soft-direction energy
curvatures are `0.006208/0.005168` for `(0.9,0)` and `0.005281/0.005594`
for `(0.9,0.3)` in `eV/(amu Å²)`; corresponding force-derived curvatures are
`0.005681/0.005685` and `0.005815/0.005819`. The largest direct
energy–force curvature discrepancy is `0.000534 eV/(amu Å²)`, well below
these sampled positive curvatures, but is not a rigorous global error bound.
The mixed direction is a metric-dependent atom–strain coordinate, **not** an
additional endpoint phonon. Its atomic and strain gradient contributions are
separately retained in the independent audit files.

All calculations preserve the five-atom `1×1×1` cubic Γ-mode source and
ABACUS PBE/`ecutwfc=100 Ry`/10 au DZP/`4×4×4` electronic k-point contract.
The latter is a k-point mesh, not a phonon supercell. No Γ re-expansion or
cutoff change is implied by this bundle.

Until a held-out interior point validates interpolation and competing branch
outcomes are mapped more broadly, these are local
conditional candidates—not a certified 2D potential-energy surface, a
T→C VCNEB barrier, or a claim of global minimum branch identity. The raw
hf workspace and its signed files remain the authoritative calculation source.
