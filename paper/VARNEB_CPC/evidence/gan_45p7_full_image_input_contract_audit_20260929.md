# GaN 45.7-GPa full-image input-contract audit (2026-09-29)

Read-only inspection of the **last visible generated inputs** in the five
completed production trees on `hf`. The exact tree paths, Slurm job IDs,
representative image-15 source hashes, and endpoint structure hashes are in
`gan_45p7_remote_source_index_20260927.md`. This audit extends the earlier
image-15 check in `gan_45p7_input_contract_audit_20260928.md`; it does not
replace the separate final-trajectory or raw-output audits.

The directory indices are contiguous: 0--28 for ABACUS, QE, ABINIT, and
CP2K; 1--27 for VASP's active interior `INCAR/KPOINTS` (fixed endpoint
statics are cached separately). VASP also retains 29 identical `POTCAR`
copies. Each count below is the number of files in which the stated setting
was observed, not a count of calculator launches or NEB iterations.

| Backend | Final visible input check | Stable potential/basis identity |
| --- | --- | --- |
| ABACUS | 29/29 `INPUT`: `basis_type lcao`, `dft_functional pbe`, `ecutwfc 100.0`; 29/29 `KPT`: Gamma 4×4×3 | 29/29 `STRU` name Ga/N UPF and the Ga/N 10-au/100-Ry DZP orbital files. For each of the four copied files, all 29 SHA-256 hashes are identical to the image-15 hashes in the source index. |
| VASP | 27/27 `INCAR`: `ENCUT 600`, `GGA PE`, `ISYM -1`, `SYMPREC 1e-4`; 27/27 `KPOINTS`: Gamma 8×8×6 | 29/29 `POTCAR` copies have SHA-256 `f94781ce6cf9b9454c093353320c5e80a2a0aceea6fb8353b33b01ac37b95168`. No POTCAR bytes are redistributed. |
| QE | 29/29 `espresso.pwi`: `ecutwfc 100`, `ecutrho 600` Ry; automatic 4×4×3, zero offset | All 29 name the same PseudoDojo NC-SR PBE v0.4 Ga/N UPF directory and filenames. |
| ABINIT | 29/29 `abinit.in`: `ecut 1400 eV`, `ngkpt 4 4 3` | All 29 name the same PseudoDojo NC-SR PBE v0.4 Ga/N PSP8 paths. |
| CP2K | 29/29 `cp2k.inp`: `CUTOFF [eV] 10884.554409746897` (= 800 Ry using the ASE Ry constant), `REL_CUTOFF 80` Ry, Monkhorst--Pack 4×4×3 | Both Ga and N use `DZVP-MOLOPT-SR-GTH` and `GTH-PBE` in all 29 files. |

The shared QE/ABINIT pseudopotential files have SHA-256 Ga/N UPF
`eb29116c4cce510e155a7e86f185d45249721df1a05ab324f9246602c8132197` /
`ae7445db618be9b28e7e273e44ad66d42eecd17cad5ab383396a5e608beb9f0f`,
and Ga/N PSP8
`8d656e2af5c22295672ec51bdefbe8c9950610fbc9387221acccf55a3d9c3259` /
`c919ba1feb97d5449a12ab14f8c97986f61e2ceb595505c46ac73628b274cb37`.
The [CP2K 2024.1 MGRID reference](https://manual.cp2k.org/cp2k-2024_1-branch/CP2K_INPUT/FORCE_EVAL/DFT/MGRID.html)
confirms that an unqualified `REL_CUTOFF` is in Ry.

For repeatability, the following SHA-256 digests commit to the lexically
ordered `sha256sum` output of the listed glob groups in each exact result
tree (`sha256sum GLOB1 GLOB2 ... | sha256sum`). The source-index tree paths
must be substituted for `TREE`. Group ordering matters; the output includes
absolute file paths. These are compact integrity commitments, not a public
archive of the calculator inputs.

| Backend | Ordered glob groups relative to `TREE` | SHA-256 of file-hash listing |
| --- | --- | --- |
| ABACUS | `*/INPUT`, `*/KPT`, `*/STRU` | `0668d71a57dcb40c2c5e6eda85797aae13dabcfb30f48ecb10a948f18f562a70` |
| VASP | `*/INCAR`, `*/KPOINTS`, `*/POTCAR` | `84cb8640f82d48207cb9b8c421455f974cb6ffcfcae6f960e0188c0384031fc9` |
| QE | `image_*/espresso.pwi` | `67d1480fd14273f0864a33af8e36f4a3e482cd6c869b9cf6ca7e4bdf799d4cc3` |
| ABINIT | `image_*/abinit.in` | `c0ffe0d38f2e8d1f50a984aaa0aac8996261cbd6d5fa594a3638cfee8bc1fe86` |
| CP2K | `image_*/cp2k.inp` | `cb14011183d2b41e284628bf198d66882b1fa4484170b8816e96dd78a67fbee9` |

This verifies static input consistency in the **last visible files**. It
does not establish that earlier overwritten calculator calls used the same
settings, that every SCF terminated normally, that potential/basis choices
are numerically converged, or that original raw energy/force/stress records
are complete. The external pressure of 45.7 GPa and endpoint-caching policy
are manager settings, not keys in these static calculator inputs. In
particular, CP2K's final-chain forces are preserved in exact-geometry caches
but original-text atomic forces and final run-end markers remain unavailable.
