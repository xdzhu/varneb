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
| B4→candidate | `27850073` | 14 | Completed; converged at FIRE step 0, final maximum generalized force 0.096528 eV/Å |
| candidate→B1 | `27850074` | 12 | Completed; converged at FIRE step 3, final maximum generalized force 0.074809 eV/Å |

Do not infer completion from Slurm disappearance alone. For each job, require a
completed manager summary, `final_max_generalized_force_eV_per_A ≤ 0.10`,
complete per-interior SCF/energy/force/stress evidence, valid geometry, and a
matching endpoint cache hash. An isolated one-step force rebound is not a
failure criterion. Recheck the stitched chain's exact common center, both
one-sided tangents, geometry, pressure convention, total enthalpy profile,
and whether any interior image overtakes the fixed center. If either segment
does not converge within wall time, restart from its latest *complete* chain
without changing the electronic contract.

For display, distinguish actual refined images inside the measured local
`(q_u,q_v)` rectangle from any interpolation between images; do not imply
that the full B4→B1 path lives inside that local chart. The full stitched
profile must use a recomputed full-path arc coordinate. Label the two fixed-
center segments and the measured DFT grid unambiguously. A visually smooth
join is not a stationary-TS certificate.

## Completed result and figure boundary

`scripts/audit_gan_45p7_two_segment_result.py` compared all 26 active-image
OUTCARs to the final complete trajectory snapshots. Each raw output contains
an electronic-convergence marker and normal VASP completion footer; energy,
force, stress, geometry, generated INCAR/KPOINTS/POTCAR hashes, and effective
symmetry contract agree with the manager record. The left and right fixed
center are the same structure and enthalpy. In the 29-image stitched chain the
center remains the highest sampled image: the barriers relative to B4 and B1
are 338.465 and 343.061 meV/GaN, respectively. These are per formula unit;
the calculation cell contains two GaN formula units.

The incoming/outgoing one-sided generalized tangent cosine at the join is
0.9636. Their squared overlaps with the local negative-curvature joint mode
are 0.9762 and 0.9863. Thus the segments meet in a similar direction, but
neither the join nor the path-force gate remedies the independently reported
normal-strain energy/stress derivative discrepancy at the candidate. Do not
call it a formally certified stationary/index-one transition state.

The measured 17×17 local mode chart spans only $|q_u|\le0.020$ Å and
$|q_v|\le0.0125$ Å. Projection of the stitched chain puts **only image 15**
inside this rectangle. Image 14 projects to approximately
$(+0.07221,+0.00944)$ Å and image 16 to $(-0.10237,+0.00965)$ Å. Therefore
the dashed center-to-neighbor rays in the new figure are clipped linear
coordinate connections, not additional NEB images, DFT statics, or a relaxed
local MEP. Panel (b) plots all 29 actual images against a recomputed full-chain
generalized arc fraction. This honest overlay is the current limit of the
fixed-center two-segment refinement. If an image-resolved line *within* the
small local mode chart is required later, near-center image insertion and
fresh same-contract DFT evaluations are necessary.

Reproducibility bundle:
`paper/VARNEB_CPC/evidence/gan_45p7_split_20261004/` (raw audit, manifest,
two trajectory histories and manager summaries, joint Hessian modes), plus
`paper/VARNEB_CPC/figures/gan_45p7_split_local_mode_20261004_v4_*` (source
table and QA record alongside PDF/SVG/PNG). Licensed POTCAR contents stay on
hf. Recreate the figure with
`python -m scripts.plot_gan_45p7_two_segment_on_local_cut` using that evidence
directory, its `joint_hessian.npz`, and the archived
`gan_600eV_local_joint_dft289_20260930_v2_source_data.csv`.
