# GaN B4→fixed candidate→B1 VASP split-chain evidence

Two non-climbing VCNEB halves were run at 45.7 GPa with the unchanged
VASP/PBE/Ga_d+N/600-eV/Γ-8×8×6 static contract. `manifest.json` records the
source-chain, source-OUTCAR and shared-center structure hashes. The complete
trajectory histories and manager summaries are archived as `left.*` and
`right.*`; the final snapshots contain 16 and 14 images, respectively, and
share one identical fixed endpoint. `raw_audit.json` records per-image hashes
and checks of all 26 interior OUTCARs on hf. Licensed POTCAR contents are not
included. `joint_hessian.npz` is the prior local Hessian-mode basis used only
for path projection and curvature diagnostics.

The raw outputs remain at the hf run root documented in
`docs/GAN_B4_B1_FIXED_CENTER_TWO_SEGMENT_2026-10-04.md`. The audit and
plotting sources are `scripts/audit_gan_45p7_two_segment_result.py` and
`scripts/plot_gan_45p7_two_segment_on_local_cut.py`. The plotted 289-node
local surface is separately archived as
`paper/VARNEB_CPC/figures/gan_600eV_local_joint_dft289_20260930_v2_source_data.csv`.

The center is the highest **sampled** chain image and a local index-one
candidate, not a certified stationary transition state. Only this center
falls inside the small measured two-mode rectangle; the plotted short dashed
rays are clipped coordinate interpolations, not newly evaluated images.
