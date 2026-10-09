# G1 preserving-band sampling audit: six registered static points

This finite audit starts from the complete ordinary-converged preserving
band of job `28319570`, step 69, archived in
`../switching_converged_update_20261009_1255/preserving_step69/`.
Its sampled maximum is 32.806023 meV/HfO2, not a certified continuous MEP
barrier. Ordinary residual 0.099372458 eV/A does not control the energy
maximum between the nine images.

Using the cached physical forces and stress, the energy derivatives along
the existing fractional-coordinate/cell line segment change from positive
to negative in segments 2->3 and 5->6. Cubic Hermite estimates suggest
38.217947 meV/HfO2 within each segment. These estimates select new samples;
they are **not additional DFT results, stationary saddles, or bounds**.

`manifest.json` freezes three fractions (0.25, 0.50, 0.75) in each segment,
six geometry hashes, source step/observation hash, physical input hashes,
common PO+ energy and the derivative evidence. Both sides are calculated,
not inferred from mirror symmetry. Structures use the original unwrapped
ordered fractional lift and linearly interpolated cell. No minimum-image
remapping, alignment or relaxation is introduced.

## Execution and cost boundary

One `hfacnormal01` allocation uses one node, 32 MPI ranks and one thread
per rank. Six SCFs run **serially** with `mpirun -np 32`, alongside at most
one other study chain. Wall-time cap is one hour. Maximum new SCF count is
six, independently recorded from the earlier G0 and chain calculations.
At the previously observed 103--154 seconds per comparable static SCF,
the provisional compute estimate is 11--16 minutes plus launcher overhead;
queue time and different geometries are not guaranteed by that estimate.

The six original INPUT/KPT/UPF/orbital hashes remain mandatory:
ABACUS/PBE, 100 Ry, complete 10-au DZP, Gamma 2x2x2, P=0 and E=0.
Only the structure changes. Ordinary NEB target stays 0.10 eV/A;
no CI, new endpoint, Hessian, chain restart or G2 submission is triggered.
Existing calculation/result directories are refused, not overwritten or
automatically repeated after failure.

Preparation is zero-DFT:

```powershell
python -m scripts.prepare_hfo2_sampling_bridge prepare --observation benchmarks/hfo2_channels/20261008/switching_converged_update_20261009_1255/preserving_step69 --root NEW_SIX_POINT_ROOT
```

The reviewed Slurm entry is
`cluster/hf_hfo2_G1_sampling_bridge_20261009.slurm`; `RUN_DFT=1`,
`SOURCE_ROOT` and `BATCH_ROOT` must be explicitly supplied. Submission
requires checking the live two-calculation capacity first.

Each point retains raw SCF output, MPI-size/convergence audit, all
energy/force/stress fields and input/structure/log hashes. The final
summary re-reads actual raw results and checks the prepared geometry and
common-reference conversion. A higher linear-reconstruction sample
indicates insufficient sampled-band resolution; it is not an optimized
continuous MEP, a relaxed conditional surface or a TS certificate.

This is an adaptive G1 sampling check, selected using already seen chain
data, **not a blind holdout or evidence for the later predictive claim**.
The existing G1->G2->G3 gates and finite material matrix are unchanged.
New submission, observed cost and audited labels will be recorded only
after they actually occur.

## Figure contract (defined before export)

Core comparison: six actual statics test whether the original sampled
maximum misses higher energies on its specified linear reconstruction.
The quantitative two-panel grid uses the existing Python/matplotlib
workflow at 183 x 99 mm. Panel (a) locates nine reused band samples and
six new SCFs on the original normalized generalized arc; panel (b)
shows the two tested segments using their dimensionless fraction.
Straight joining lines are sample guides, not a smooth MEP fit. A dotted
reference marks the original sampled maximum, not an error bound.

Both panels share the energy scale, aligned axes and normal-weight
larger `(a)`/`(b)` labels. All four spines are visible; framed translucent
white legends remain outside curves. Export editable SVG/PDF, 600-dpi
TIFF, PNG preview, all 15 source rows and numeric/provenance QA. Review
actual pixels before delivery. There are no independent replicates,
statistical error bars or stationary-TS labels. The review risk is
mistaking the higher **linear reconstruction** samples for the peak of
an independently relaxed continuous MEP, or calling this adaptive audit
a blind prediction. Those interpretations are explicitly excluded.

## Actual result: job 28374431, completed 2026-10-09

Slurm reports COMPLETED0:0, 13:52:25--14:02:26 CST, one node/32 CPUs.
All six raw SCFs independently pass convergence, DSIZE32, complete
finite energy/forces/stress, input-hash and STRU round-trip checks on HF.

| Segment | Fraction | Energy above common PO+ (meV/HfO2) |
|---|---:|---:|
| 2->3 | 0.25 | 34.4379586541 |
| 2->3 | 0.50 | 38.5084306490 |
| 2->3 | 0.75 | 37.9067777562 |
| 5->6 | 0.25 | 37.9067616818 |
| 5->6 | 0.50 | 38.5084480845 |
| 5->6 | 0.75 | 34.4379828925 |

Highest new sample38.5084480845 exceeds the original sampled maximum by
5.7024248676meV/HfO2. It is a reconstruction sample, not a certified TS.
The two sampled maxima differ by0.0000174355meV/HfO2; this mirror agreement
does not estimate full barrier error or functional uncertainty.

The original nine cached images plus the six real statics form a fifteen-image
linear reconstruction (thirteen moving images). Its ordinary residual
replays0.09937245841139401eV/A: **no optimizer steps, no additional SCFs**.
There is no automatic new chain submission. This shows directly why an
ordinary force pass alone cannot resolve the energy maximum between images.

SCF elapsed times sum542.979671s (74.88--99.36s per SCF), or4.826486core-hours
at32ranks. Allocation wall time601s gives5.342222core-hours including overhead.
`completed_HF/` preserves the actual manifest, summary, rank probe, geometry,
E/F/stress, STRU, INPUT/KPT, physical-input hash records and raw logs. Repeated
large UPF/orbital payloads are not copied into the repository; their exact
bytes were checked on HF, not reconstructed from the records.

The local ASE distribution has no ABACUS I/O plugin. Offline geometry review
therefore replays the **same fixed production STRU writer** and compares every
token (nonnumeric tokens exactly; floating text tails within1e-13), in addition
to raw STRU/geometry hashes. This is separate from the actual HF ASE-ABACUS
round-trip checks, not an installation or change to the production calculator.

The new [two-panel figure](../../../../paper/VARNEB_JCTC/figures/hfo2_G1_sampling_bridge_20261009/hfo2_G1_sampling_bridge.png)
replays all source energies and raw logs before drawing. Reproduce from
the repository root:

```powershell
python -m scripts.plot_hfo2_sampling_bridge --observation benchmarks/hfo2_channels/20261008/switching_converged_update_20261009_1255/preserving_step69 --completed benchmarks/hfo2_channels/20261008/preserving_sampling_bridge_20261009/completed_HF --output NEW_FIGURE_BUNDLE
```

Both the local preparation and independently generated HF preparation remain
archived: their cell matrices are identical and unwrapped fractional
coordinates differ at most2.22e-16. Exact physical template hashes remain
unchanged; these generated geometry floating tails are not electronic-input
retuning. See `preparation_validation.json` and the final delivery receipt.
