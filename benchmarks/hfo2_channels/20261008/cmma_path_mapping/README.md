# Cmma directions registered to the production ordered-path chart

2026-10-09, E023. Zero new DFT. This extends the
[published-reference internal audit](../cmma_reference/README.md) to existing
PBE/ABACUS paths. It is a descriptive strong-reference preparation, not a
successful barrier predictor, conditional surface or TS certification.

## Registration must precede coverage evaluation

The published Cmma and production T cells cannot be identified by their axis
labels alone. The author's T distortion axis differs from our T axis. A direct
nearest-Cmma-to-T mapping also contains four equal-distance oxygen assignments.
The initial arbitrary-axis trial was therefore rejected as a strong-control
comparison, not used to conclude that the reference fails.

The fixed procedure first enumerates 24 proper signed axis frames and eight
quarter-cell origins (192 candidates). Register the author's pinned T structure
against our stationary T in the original T metric, keeping all cost ties within
`1e-7 Å`; 16 frames survive. Among these, register Cmma to the **known common
PO+ endpoint**, retaining all four tied minima. Each selected species assignment
has a global second-best cost gap of `1.47265 Å`. The complete geometric scores,
integer shifts, proper frames and assignments are retained in [analysis.json](analysis.json).
No energy, mode coverage, force or strained holdout enters this choice.
This is G1-informed exploration, not a blind or uniquely physical atom mapping.

The additional frame anchor is the author's
[Tetragonal.vasp](https://raw.githubusercontent.com/yuboqiuab/unstableflatband/a438e4ecf63cddfa68d4ae8ea83d4c19bd238969/figure1/Tetragonal.vasp),
919 bytes, Git blob `b8cb27541fa7d1fe166b37fa9c432d31cf664047`, SHA256
`afcb23ec69f99d4ef9fa546795222d4c28eac4249659eb35c528fe5f82131a53`.
Its pinned public GitHub metadata was checked against the original file.

Every path image keeps its original order, cell, coordinates and continuous
periodic lift. The selected endpoint assignment maps **reference rows only**;
no image-wise rematching, MIC folding, structure repair or optimizer change occurs.
All equivalent registered frames remain in the report rather than selecting one
by favourable coverage.

## A common chart and an explicit real-span policy

With ASE row cells, rotated reference displacements enter the original production
T fractional-displacement metric through `v_common = v_Cmma @ inv(H_Cmma) @ H_T`.
This declared coordinate transport is not a physical strain or SCF input.
Masses remain Hf178.49/O15.9994 amu. Cell strain is separate from atomic coordinates.

The full four resolved negative conventional-Gamma directions are retained.
`vcneb.mode_subspaces.real_mode_subspace` performs an explicit SVD of
`[Re(V), Im(V)]`, requires declared real rank four and records complex-column
reconstruction error. The rank/projection thresholds are `5e-5`; actual maximum
relative projection error is `2.084e-6`. The real basis is orthonormal, but its
SVD axes are **not new phonon eigenpairs or named individual modes**. Subspace
weights, not arbitrary SVD coordinates, are compared. Complex-unitary changes
of a real span preserve its projector in synthetic tests.

Source rigid-translation admixture is measurable: the common-chart subspace
overlap Frobenius norm is `0.0036854`. It is retained and disclosed, not removed
by an ASR operation or dismissed as the printed Gram defect. Only displacement
analysis removes an overall mass-weighted rigid translation, as in the existing
projection contract. Raw/transformed author eigenvector arrays are not redistributed.

## What the first converged band actually establishes

Ten frozen image records of ordinary T→PO are checked against ordered E/F/stress
and geometric hashes and replayed at `fmax=0.059881615690252875 eV/Å`. In the
**same T-origin, same mass-metric displacement chart**:

| Representation | Rank | Fraction at highest image3 | Fraction at PO |
|---|---:|---:|---:|
| Cmma four-direction span, all registered frames | 4 | 0.271810–0.271901 | 0.522401–0.522402 |
| Rotated T geometric triplet | 3 | 0.896905 | 0.304151 |
| Lowest two complete T optical doublets | 4 | 0.035448 | recorded in JSON |

At zero T displacement the fraction is undefined; the API zero flag is retained,
not interpreted as measured zero coverage. The T-triplet rank differs from the
full Cmma rank. These numbers are not the older fluorite-origin/unit-metric
fractions and must not be mixed with those denominators.

Separately, the Cmma-origin affine-plane residual in this same T metric is
`6.43182 sqrt(amu) Å` at T and `5.60197 sqrt(amu) Å` at PO. Thus this *linear
four-direction plane* does not contain both endpoints. Its origin-specific
fractions have different denominators and are labelled separately. This does
not rule out nonlinear relaxation, larger fixed spaces, other conventional-cell
periodicities or another effective fixed-parent model.

[historical_observations.json](historical_observations.json) applies the **same
registration and basis policy** to nine earlier complete snapshots: 84 frozen
image records, including reused endpoints. This is not 84 new or necessarily
unique SCFs, and none of those snapshots is a newly converged competing path.
At their recorded maxima, the four Cmma directions capture about44.5% versus
46.6% for the T-triplet in PO→M step10,18.1% versus46.5% in the preserving
flip step10, and44.8% versus29.2% in the reversing flip step10. These remain
unconverged, metric-dependent descriptive observations, not a final channel
ranking or demonstration of improved predictions. They motivate evaluating
representation across channels, not only the T→PO formation path.

## Reproduction and remaining gates

Keep the six pinned author files outside Git; no author executable is run.
Use a **new output path**:

```text
python -m scripts.analyze_hfo2_cmma_path_reference --source-root AUTHOR_FILES --dataset benchmarks/hfo2_channels/20261008/converged_gap --variants benchmarks/hfo2_channels/20261008/reference_variants --gamma benchmarks/hfo2_channels/20261008/gamma_analysis/T_d0.01.npz --output NEW_REPORT.json
```

For historical snapshots, change only `--dataset` to the existing
`benchmarks/hfo2_channels/20261008/chain_observations`. No calculator is called.
See [validation_delivery.json](validation_delivery.json) for clean-tree and
same-archive HF checks and actual production-job state.

Production-gauge registration is now implemented for these frozen observations.
This does not complete G1's competing paths or G2/G3's matched-boundary material
and independent-prediction gates. Strong frozen/relaxed/branch-aware models must
still use matched training data and costs; apparent geometric compactness is
not a barrier-accuracy claim. Source LDA curvatures, energies and eigenpairs are
not imported as PBE response data. No cutoff, orbital, k mesh, SCF setting,
pressure, new phonon matrix, chain or CI change has been made.
