# GaN B4→candidate→B1 fixed-center VCNEB refinement

This is a bounded follow-up to the converged 29-image, 45.7-GPa VASP/PBE GaN
B4→B1 chain. Its shared center is the independently audited one-step-refined
image-15 static, **not** a certified stationary transition state. Splitting
the path cannot itself cure the known volume-strain energy/stress derivative
discrepancy or prove an index-one saddle.

## Input and execution contract

- Source chain: `paper/VARNEB_CPC/evidence/gan_45p7_final_chains_20260927/gan_vasp_45p7_final_chain.traj`, SHA-256 `952298c1830b293690b8fc4722649147cba236e50e149d07dae45ff75fc2bb7e`.
- Shared center: completed static job `27792422`, OUTCAR SHA-256 `66e9ddf2bf5ee1f3f510bf71b8ec8dea2e05ae3d8be60e9593e741a87dabcfe5`.
- Original B4/B1 static OUTCAR SHA-256: `3c6795acf28dad6d00dd33a6f99cd755bf0b0cf8c9c4490bcc1fb3de3005dc87` / `5b69743b754b170807ada26d9ed147de3deb83fc6bc7b8749434b2f7cf98a671`.
- Electronic settings are unchanged: VASP 6.3.2, PBE Ga_d+N PAW, ENCUT 600 eV, Γ-centered 8×8×6, EDIFF 1e-7, ISYM=-1, SYMPREC=1e-4. External pressure is 45.7 GPa. The licensed POTCAR remains on hf.
- Left segment retains source images 0–15: 16 total, 14 active interiors. Right retains 15–28: 14 total, 12 active interiors. Only their common fixed endpoint replaces the archived chain's image 15. All three fixed structures reuse completed raw statics through geometry- and input-checked endpoint caches; they are not recalculated on each VCNEB step.
- Both runs use ordinary, non-climbing VCNEB with the source path's MIC convention, FIRE, spring `k=0.20`, and generalized-force threshold 0.10 eV/Å. Each manager uses 2 hfacnormal01 nodes and six concurrent 32-MPI image workers (194 reserved tasks including spare manager capacity), with a 96-hour wall limit.

The preparation script is `scripts/prepare_gan_45p7_two_segment_vcneb.py`; the job template is `cluster/hf_gan_45p7_two_segment_vcneb.slurm`. Preparation and both VASP `--validate-only` preflights passed without DFT. The archive/execution root is:

`/public/home/iai806/abacus/agent-runs/20261004-varneb-gan-split-600eV`

The prepared case is `case-v2`. An initial `case` preparation stopped before its manifest was written because NumPy calculator parameters were not JSON-serializable; `case-v2` includes the fixed serializer and was fully preflighted. No DFT was launched from the incomplete preparation.

## Submitted jobs and acceptance checks

| Segment | Slurm job | Expected active interiors | Submission state |
|---|---:|---:|---|
| B4→candidate | `27850073` | 14 | Running at initial check, 2026-10-04 |
| candidate→B1 | `27850074` | 12 | Running at initial check, 2026-10-04 |

Do not infer completion from Slurm disappearance alone. For each job, require a
completed manager summary, `final_max_generalized_force_eV_per_A ≤ 0.10`,
complete per-interior SCF/energy/force/stress evidence, valid geometry, and a
matching endpoint cache hash. An isolated one-step force rebound is not a
failure criterion. Recheck the stitched chain's exact common center, both
one-sided tangents, geometry, pressure convention, total enthalpy profile,
and whether any interior image overtakes the fixed center. If either segment
does not converge within wall time, restart from its latest *complete* chain
without changing the electronic contract.

For display, project only the refined images inside the measured local
`(q_u,q_v)` rectangle onto the 289-point cut; do not imply that the full
B4→B1 path lives inside that local chart. The central transverse `(s,q⊥)`
figure can show the stitched path with an explicitly recomputed full-path arc
coordinate. The figure legend must distinguish the original production chain,
the two fixed-center refined segments, and the measured DFT grid. A visually
smooth join is not a stationary-TS certificate.
