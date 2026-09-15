# VARNEB claim-to-evidence map

This map keeps the manuscript claims tied to executable code, regression tests,
and archived calculation evidence. A material result is only treated as
production evidence when its endpoint policy, calculator settings, worker
manifest, trajectory summary, and independent audit are available.

| Manuscript claim | Implementation / test | Material or model evidence |
| --- | --- | --- |
| Fractional atomic coordinates and deformation-gradient cell force are finite-difference consistent | vcneb/core.py; tests/check_vcneb_forces.py (force and pressure variants) | Analytic coupled atomic--cell model in the same regression |
| Fixed-cell reduction and ASE-compatible optimizer interface | vcneb/core.py; tests/check_vcneb_forces.py (fixed-cell comparison) | outputs/vcneb_p0_baseline_manifest.json |
| Mode-guided, strict/projected mode, and direction-basis semantics are explicit | vcneb/modes.py, vcneb/core.py; tests/check_vcneb_forces.py mode and direction regressions | examples/compare_mode_path_variants.py; docs/theory.md |
| Calculator contract requires finite energy, forces, and stress before VC-NEB | vcneb/calculator.py; tests/check_vcneb_forces.py calculator-contract regression | outputs/abacus_hfo2_smoke_manifest.json; outputs/vasp_gan_single_image_manifest.json |
| Image workers run concurrently while the manager owns the path | vcneb/executor.py; tests/check_vcneb_forces.py threaded executor and endpoint-ownership regressions | HfO2 worker manifests record interior indices only and 32-MPI workers |
| Exact-state caching and restart preserve successful sibling work after a failed batch | vcneb/executor.py; tests/check_vcneb_forces.py cache/recovery regressions | Per-run image_worker_manifest.jsonl and snapshot directories |
| Failure reports are atomic and classify recoverable failure classes conservatively | vcneb/core.py, vcneb/calculator.py; tests/check_vcneb_forces.py failure-report/classification regressions | vcneb_failure.json emitted by the ABACUS/VASP examples on interruption |
| HfO2 ordinary 7/9-image paths provide an image-count check | scripts/compare_vcneb_images.py, scripts/export_vcneb_metrics.py; tests/check_image_comparison.py, tests/check_vcneb_metrics.py | outputs/hfo2_t_to_po_pbe100_dzp10au/hfo2_validation_provenance.json, hfo2_vcneb_closure.md, image_count_comparison.md, ordinary_image_comparison.json; jobs 27678218 and 27678407 |
| HfO2 cell interpolation sensitivity is separated from image-count convergence | scripts/compare_vcneb_images.py (`--allow-duplicate-image-counts`); tests/check_image_comparison.py | outputs/hfo2_t_to_po_pbe100_dzp10au/cell_interpolation_comparison.json; job 27678924 and archived linear-cell manifest |
| Key structural changes and modal coordinates can be inspected without rerunning DFT | scripts/export_vcneb_structural_metrics.py; tests/check_geometry_metrics.py | outputs/hfo2_t_to_po_pbe100_dzp10au/vcneb_n7_fire_distributed_linear_job27678924/vcneb_structural_metrics.csv; endpoint-displacement projection is explicitly not a phonon eigenvector |
| HfO2 CI refinement is only run after a real ordinary interior peak and survives a rebound window | run_vcneb staged-CI semantics; docs/theory.md | outputs/hfo2_t_to_po_pbe100_dzp10au/vcneb_n7_ci_refine_job27678507/; job 27678507 |
| BaTiO3 T-to-C is barrierless for the tested image counts | tests/check_vcneb_audit.py; scripts/export_vcneb_metrics.py | outputs/batio3_t_to_c_pbe100_dzp10au/bto_validation_provenance.json, bto_convergence_matrix.md; CI gate is explicitly withheld |
| Public package metadata and release policy are deliberate | pyproject.toml, .github/workflows/publish-pypi.yml; tests/check_release_metadata.py | Isolated python -m build produces varneb-0.0.1; the next description update requires a new version because PyPI releases are immutable |

## Claim boundaries

The HfO2 and BaTiO3 production evidence is an ABACUS/PBE family, not a
universality proof. VASP is currently validated at single-image smoke-test
level only. A five-total-image HfO2 control was deliberately retained as
non-converged and is not used to quote a barrier. The independent linear-cell
interpolation job 27678924 is a completed, audited sensitivity comparison and
is not folded into the image-count spread. The endpoint-displacement
mode-guided HfO2 branch through job 27687189 is geometry-valid but still above
the force threshold, so its small apparent barrier is retained only as an
unconverged diagnostic.
