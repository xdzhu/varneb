# GaN local joint-mode surface: 9×9 → 17×17 refinement

The manuscript's existing 81-point local atom–strain cut is sampled on a
9×9 rectangle, `q_u ∈ [−0.020, 0.020] Å` and
`q_v ∈ [−0.0125, 0.0125] Å`. The earlier 5×5 interpolant predicted the
56 new 9×9 DFT nodes with a 0.01349-meV/GaN maximum absolute error. This
is already a small *interpolation* error, but the user requested a visibly
denser grid of genuinely calculated points. The natural nested step is a
17×17 grid: all 81 existing nodes are reused and 208 midpoints are new.
No enlarged coordinate domain or parameter change is implied.

## Immutable physical contract

- VASP 6.3.2, PBE, Ga_d+N PAW, 600-eV ENCUT, Γ-centered 8×8×6 mesh,
  `ISYM=-1`, `SYMPREC=1e-4`; the same four input-file hashes as the 9×9 cut.
- Enthalpy uses `H=E+PV` at 45.7 GPa, normalized by two GaN formula units
  per four-atom cell. The center, Hessian eigendirections, coordinate signs,
  cell/atomic displacement convention, and rectangle are unchanged.
- Every new node is one frozen-geometry VASP static. It is **not** a
  conditional relaxation. This remains a local cut, not a full-path 2D
  surface, a lower MEP, or an index-one transition-state certificate.

## Sampling and acceptance

`scripts/prepare_gan_600eV_ts_2d_dense17.py` verifies the original pilot,
5×5, 9×9, Hessian, and center hashes before writing 208 inputs. It rejects
duplicate/off-grid coordinates, Bravais-risk geometry, failed POSCAR
roundtrips, or any electronic-input hash mismatch. The prospective test is
declared **before** seeing the new DFT energies: interpolate only the old
81-point surface and score its predictions at all 208 new coordinates.
For a polished contour, require maximum absolute error ≤0.02 meV/GaN and
RMS error ≤0.01 meV/GaN. The audited 17×17 contour will still display all
measured nodes and explicitly label interpolated pixels as visual guides.
If either gate fails, report the failure and inspect its location rather
than hide points or silently adjust the calculation.

The 56 tasks in the previous 9×9 Slurm array took 118–156 s each, with a
135-s median using 32 MPI ranks. At four concurrent tasks, 208 analogous
statics represent roughly 1.9 h of pure execution and about 250 core-hours;
queue waits and outliers are additional. The account's simultaneous-submit
limit is 200 shared with other work, so the 208 points are divided into
disjoint Slurm arrays. `27807980` covers `0–103`; `27808164` covers
`104–155`; `27808203` covers `156–169`. The final `170–207` indices
remain **unsubmitted** until quota is available. All arrays use the same hash-pinned staged manifest and
`cluster/hf_gan_600eV_ts_2d_dense17.slurm`. The scheduler's test-only
start dates have been unreliable; use actual `squeue`/`sacct` state for
progress and completion estimates.

After both arrays finish, run the raw audit
`scripts/audit_gan_600eV_ts_2d_dense17.py` and only then render the new
figure with `scripts/plot_gan_600eV_local_joint_cut17.py`. The plot's right
panel will compare the 208 genuinely new statics with the frozen prior
9×9 predictions, while the left panel overlays all 289 measured nodes.
Only an audited, gate-qualified result may replace the current 81-point
Fig. 5 in the CPC manuscript.
