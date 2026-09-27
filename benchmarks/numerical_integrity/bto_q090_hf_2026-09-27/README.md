# BTO off-axis conditional-branch evidence, `Q_z=0.9`

These sixteen JSON files are byte-for-byte copies of the completed `27786180`
and `27787086` ABACUS runs and subsequent read-only audits on hf. Their source directory is
`/public/home/iai806/abacus/agent-runs/20260927-varneb-bto-soft-gridpilot-qx030`.
The `run-q090_q000/` and `run-q090_q030/` subdirectories there retain the
original `INPUT`, `KPT`, `STRU`, ABACUS logs, and immutable evaluation cache.
The local copies here retain the SHA-256 links back to those raw inputs and
logs; they do not duplicate the large raw calculation tree.

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

All calculations preserve the five-atom `1×1×1` cubic Γ-mode source and
ABACUS PBE/`ecutwfc=100 Ry`/10 au DZP/`4×4×4` electronic k-point contract.
The latter is a k-point mesh, not a phonon supercell. No Γ re-expansion or
cutoff change is implied by this bundle.

Until the selected branches pass independent curvature audits at a second step,
direct soft-direction energy/force checks, and a held-out interior point
validates interpolation, these are local
conditional candidates—not a certified 2D potential-energy surface, a
T→C VCNEB barrier, or a claim of global minimum branch identity. The raw
hf workspace and its signed files remain the authoritative calculation source.
