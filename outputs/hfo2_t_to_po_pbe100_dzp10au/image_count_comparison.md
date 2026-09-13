# HfO2 T→PO image-count comparison

All production branches use ABACUS PBE, 100 Ry, full 10 au Orb-DZP basis sets for the elements in each material, 2×2×2 k-points, `k=0.2`, FIRE, and fixed endpoint energies/forces/stress cached once. `N-image` means the total number of frames, including the two endpoints; only `N-2` interior frames are dispatched to 32-MPI workers.

| branch | Slurm job | total / interior | status | final max generalized force (eV/Å) | forward enthalpy barrier (eV) | reaction enthalpy (eV) | peak image |
|---|---:|---:|---|---:|---:|---:|---:|
| ordinary | 27678406 + 27678689 | 5 / 3 | not converged; cancelled after extended rebound window | 0.098970 best retained (0.112725 in controlled continuation) | not claimed | not claimed | — |
| ordinary | 27678218 | 7 / 5 | completed | 0.0455156 | 0.1567509 | -0.3252848 | 2 |
| ordinary | 27678407 | 9 / 7 | completed | 0.0497191 | 0.1562092 | -0.3252848 | 3 |
| CI refinement | 27678507 | 7 / 5 | completed after rebound/platform observation | 0.0295475 | 0.1291722 | -0.3252848 | 2 (climbing) |

The converged ordinary 7/9-image barriers differ by only `0.0005417 eV` (about 0.35% relative to the 7-image value), and the reaction enthalpy is identical to the displayed precision. The shift of the discrete peak from image 2 to 3 is expected when the path is sampled more finely; it is not evidence for a second saddle. The 5-image branch is excluded from the physical barrier comparison because its generalized force remained above target and continued to rebound even after a smaller-step continuation.

Each completed 7/9-image branch passed `scripts/audit_vcneb_result.py --max-min-distance 2.0` with `status=ok` and no issues. The copied machine-readable summaries, audits, preflight reports, and worker manifests are in:

- `vcneb_n7_ci_refine_job27678507/`
- `vcneb_n9_fire_distributed_cmp_job27678407/`

The authoritative full ABACUS scratch directories remain on `hfacnormal01` under `/public/home/iai806/abacus/agent-runs/20260912-vcneb-hf/validation/hfo2_t_to_po/`.
