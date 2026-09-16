# Material validation and figure-reproduction guide

This guide connects the two completed material cases to the calculator-free
reporting tools. It does not launch a calculator, submit Slurm work, or alter
an archived trajectory.

## What each case establishes

| Case | Scientific role | Accepted conclusion | Not established |
| --- | --- | --- | --- |
| BaTiO3 T to C | compact, direct perovskite path | 5/7/9 total images are monotonic and consistently barrierless; CI is withheld | a finite-temperature transition rate or a first-order saddle |
| HfO2 T to PO | 12-atom reconstructive, variable-cell path | ordinary image-count stability, linear/log-strain control, and a staged CI barrier | universality over mappings, supercells, strain states, or calculators |

All production results use ABACUS/PBE and must be described as that calculator
family. The VASP adapter has a real single-image smoke test but no production
material path, so it cannot be presented as a cross-engine barrier benchmark.

## Regenerate the paper figure from archived summaries

Run from the repository root after the completed summary JSON files are
available locally:

```text
python scripts/plot_material_validation_figure.py \
  --bto "5 images=outputs/batio3_t_to_c_pbe100_dzp10au/vcneb_t_to_c_n5_parallel_converged_job27676251/vcneb_summary.json" \
  --bto "7 images=outputs/batio3_t_to_c_pbe100_dzp10au/vcneb_t_to_c_n7_parallel_converged_job27676310/vcneb_summary.json" \
  --bto "9 images=outputs/batio3_t_to_c_pbe100_dzp10au/vcneb_t_to_c_n9_parallel_converged_job27676513/vcneb_summary.json" \
  --hfo2 "ordinary n7 log-strain=outputs/hfo2_t_to_po_pbe100_dzp10au/vcneb_n7_fire_distributed_job27678218/vcneb_summary.json" \
  --hfo2 "ordinary n9 log-strain=outputs/hfo2_t_to_po_pbe100_dzp10au/vcneb_n9_fire_distributed_cmp_job27678407/vcneb_summary.json" \
  --hfo2 "CI n7=outputs/hfo2_t_to_po_pbe100_dzp10au/vcneb_n7_ci_refine_job27678507/vcneb_summary.json" \
  --hfo2 "ordinary n7 linear-cell=outputs/hfo2_t_to_po_pbe100_dzp10au/vcneb_n7_fire_distributed_linear_job27678924/vcneb_summary.json" \
  --output-dir outputs/vcneb_material_validation_figure
```

The command writes a PNG preview, editable SVG/PDF, and two CSV source-data
tables. It rejects incomplete summaries and requires an explicitly labeled
CI path. The figure normalizes each path by its archived geometric arc length;
it does not imply that images with different counts sit at identical physical
configurations.

## Reporting rules

1. State whether image count means total images or inserted/interior images.
   Here, 7 total images means 5 interior images and fixed cached endpoints.
2. For BTO, report a monotonic endpoint rise, not an internal activation
   barrier. A CI calculation is inappropriate without an interior peak.
3. For HfO2, keep ordinary 7/9 image-count comparison, linear-cell control,
   and CI refinement as distinct rows. Do not fold the linear-cell control
   into the image-count spread.
4. State force semantics. New ordinary NEB runs use the default maximum
   generalized-force target of `0.10 eV/A`. The historical strict production
   table uses recorded `0.05 eV/A` (and CI `0.03 eV/A`) targets, which must not
   be relabeled as current defaults.
5. Quote HfO2 values per 12-atom cell or per formula unit, but never change
   denominator inside one comparison. Four HfO2 formula units occupy the
   present cell.
6. Label the LDA/QE/USPEX 40-image HfO2 literature value as external
   comparison data, not a replica or error bar. The BTO PBEsol restrained
   NEB result is an energy-scale comparison only.

## Checks before release

```text
python tests/check_material_validation_figure.py
python -m pytest -q
```

Inspect the SVG after generation to confirm text remains editable, then keep
the source-data CSV files beside the figure in any manuscript archive.
