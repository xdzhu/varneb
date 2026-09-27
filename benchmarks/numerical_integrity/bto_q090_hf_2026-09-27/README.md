# BTO off-axis conditional-branch evidence, `Q_z=0.9`

These ten JSON files are byte-for-byte copies of the completed `27786180`
ABACUS runs and subsequent read-only audits on hf. Their source directory is
`/public/home/iai806/abacus/agent-runs/20260927-varneb-bto-soft-gridpilot-qx030`.
The `run-q090_q000/` and `run-q090_q030/` subdirectories there retain the
original `INPUT`, `KPT`, `STRU`, ABACUS logs, and immutable evaluation cache.
The local copies here retain the SHA-256 links back to those raw inputs and
logs; they do not duplicate the large raw calculation tree.

The two selected `-Q_y` branches satisfy the stated orthogonal-gradient and
1-kbar residual-stress gates. Independent raw-output audits cover 77 and 78
DFT evaluations. Three-start cache replay resolves both the higher-energy
frozen branch and the nearly degenerate `+Q_y`/`-Q_y` pair. The `0.05`
`sqrt(amu) Å` curvature-preflight files certify only that all 32 signed
probe geometries are safe; they are **not** curvature results.

All calculations preserve the five-atom `1×1×1` cubic Γ-mode source and
ABACUS PBE/`ecutwfc=100 Ry`/10 au DZP/`4×4×4` electronic k-point contract.
The latter is a k-point mesh, not a phonon supercell. No Γ re-expansion or
cutoff change is implied by this bundle.

Until the selected branches pass independent curvature audits at two steps
and a held-out interior point validates interpolation, these are local
conditional candidates—not a certified 2D potential-energy surface, a
T→C VCNEB barrier, or a claim of global minimum branch identity. The raw
hf workspace and its signed files remain the authoritative calculation source.
