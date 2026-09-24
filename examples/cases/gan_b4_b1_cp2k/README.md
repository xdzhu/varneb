# GaN B4→B1 with CP2K at 45.7 GPa

This is the compact input fixture for the independently converged GaN
tetragonal-route case. `B4_CONTCAR.vasp` and `B1_CONTCAR.vasp` are four-atom,
two-formula-unit endpoints with identical Ga/N order. They are **not** a
license to skip fresh endpoint force/stress checks on a different machine or
CP2K build.

| Contract | Production value |
| --- | --- |
| XC and potentials | PBE, GTH-PBE |
| Basis | DZVP-MOLOPT-SR-GTH |
| CP2K cutoff / relative cutoff | 800 / 80 Ry |
| Brillouin-zone mesh | 4×4×3 |
| Pressure | 45.7 GPa hydrostatic |
| Band | 29 total images, 27 interior, no climbing image |
| Spring / force criterion | 0.20 / 0.10 eV/Å |
| Mapping/gauge | identity, no cell or translation alignment, linear cell interpolation |
| Accepted result | Slurm `27770714`, step 52, `fmax=0.09409 eV/Å`, barrier `0.29276 eV/GaN`, peak image 15 |

The calculator settings are in
[`../../material_profiles/cp2k_gan_pbe_dzvp_k4_cutoff800_rel80.json`](../../material_profiles/cp2k_gan_pbe_dzvp_k4_cutoff800_rel80.json).
The evaluated final band and exact normalization audit are kept in the
manuscript's `evidence/` directory; the path is plotted in its five-backend
figure and as a separate final-path preview.

## Run sequence

1. Load the site-provided CP2K 2024.1/Intel MPI environment inside a Slurm
   allocation. Verify that CP2K resolves the GTH basis/potential data.
2. Run `examples/run_vcneb_ase.py` with `--static-only` for both endpoints
   under the fixed profile. Check the resulting `ase_static_summary.json` with
   `scripts/validate_ase_static_gate.py --pressure-gpa 45.7 --fmax 0.10
   --stress-kbar 1.0`. The accepted endpoint residuals were 0.416 and
   0.249 kbar; a new installation must establish its own passing gate.
3. Run the same driver with the matching `--endpoint-static-summary`,
   `--n-images 29 --pressure-gpa 45.7 --k 0.20 --fmax 0.10 --no-climb
   --mapping identity --no-align-cells --cell-interpolation linear
   --maximum-cell-step 0.05 --candidate-step-retries 8`, and give every image
   its own working directory. First use `--validate-only` to check the actual
   endpoint hashes and initial geometry without launching DFT.
4. On completion, require `vcneb_summary.json` to report `converged`, check
   every image energy/force/stress and geometry, and recompute
   `H_i=E_i+45.7 GPa×V_i` divided by two GaN units. Do not infer convergence
   solely from a zero Slurm exit code.

### MPI isolation warning

`IMAGE_WORKERS × CP2K_MPI_RANKS ≤ allocated tasks` is necessary but **not
sufficient**. Five independent `mpirun -np 16` commands in the same Slurm
allocation can all pin to the same first 16 CPUs. That happened in an
independent step-33 continuation (`27775731`), making each electronic step
abnormally slow; that branch was cancelled and is not an acceleration
benchmark. Use one worker until a site-specific five-worker affinity canary
demonstrates disjoint CPU sets, or provide an explicitly verified launcher.
Changing the optimizer cannot repair an MPI placement error.

The literature guide for this **tetragonal** route is approximately
`0.34 eV/GaN`; the `0.39 eV/GaN` guide applies to the separate hexagonal
route. Differences between codes are not a total-energy identity test because
basis sets, pseudopotentials, and relaxed endpoints differ.
