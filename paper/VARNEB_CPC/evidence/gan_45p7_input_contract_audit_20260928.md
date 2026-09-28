# GaN 45.7-GPa production input-contract check (2026-09-28)

This is a **representative generated-input** check for Table 1 of the CPC
manuscript, not a convergence study or a full 29-image raw-output audit. All
five samples are the actual image-15 calculator inputs from the completed
production trees indexed in `gan_45p7_remote_source_index_20260927.md`.
The ABACUS and VASP files were inspected read-only on `hf`; their remote
SHA-256 hashes were rechecked on 2026-09-28 against that index. The QE,
ABINIT, and CP2K files are byte-preserved locally in
`gan_45p7_final_peak_inputs_20260927/`; their local hashes match the index.
No VASP POTCAR is redistributed.

| Backend | Actual input evidence | Observed cutoff | Observed k mesh |
| --- | --- | --- | --- |
| ABACUS | image-15 `INPUT` (`dc6684ff…`) and `KPT` (`bda588d9…`) | `ecutwfc 100.0` Ry with 10-au DZP orbitals identified by the source-indexed `STRU` and orbital hashes | `Gamma`, `4 4 3` |
| VASP | image-15 `INCAR` (`83ba34d4…`) and `KPOINTS` (`b5215f3e…`) | `ENCUT = 600` eV | `Gamma`, `8 8 6` |
| QE | archived `espresso.pwi` (`b558fff1…`) | `ecutwfc = 100`, `ecutrho = 600` Ry | automatic `4 4 3`, zero offset |
| ABINIT | archived `abinit.in` (`b4d43889…`) | `ecut 1400 eV` | `ngkpt 4 4 3` |
| CP2K | archived `cp2k.inp` (`935170d7…`) | `CUTOFF [eV] 10884.554409746897` = 800 Ry using the ASE Ry constant; `REL_CUTOFF 80` Ry | `MONKHORST-PACK 4 4 3` |

The [CP2K 2024.1 MGRID reference](https://manual.cp2k.org/cp2k-2024_1-branch/CP2K_INPUT/FORCE_EVAL/DFT/MGRID.html)
specifies Ry as the default unit of an unqualified `REL_CUTOFF`.

The same QE/ABINIT/CP2K input files explicitly select the PseudoDojo
NC-SR PBE pseudopotential paths or CP2K GTH-PBE and DZVP-MOLOPT-SR-GTH.
The remote index pins the VASP POTCAR and ABACUS orbital/pseudopotential
hashes without shipping them. The 45.7-GPa pressure and 29-image policy
belong to the **manager** manifest, not to these single-point calculator
inputs; see `gan_45p7_multibackend_vcneb_20260924.json`.

This audit does not prove that all generated image inputs are identical in
settings, that all raw calculator calls ended normally, or that basis and
k-mesh errors are converged. The final-chain geometry/energy/force/stress
audit is separately recorded in `gan_45p7_final_chain_audit_20260927.json`.
For CP2K, the exact-geometry cache contains final forces and all 27
interior raw-text energies and stresses reconcile, but original-text atomic
forces and final run-end markers are unavailable. Those evidence tiers must
not be merged into a stronger raw-log certification.
