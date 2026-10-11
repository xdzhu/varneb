# E077: external work from the existing clamped M band

Zero DFT: reparse the nine portable native E068 INPUT/KPT/STRU/log/E/F/stress
records. The original HF audit, not this workstation, attests the full six
physical files and91fresh SCFs. The ordinary residual remains0.099214088.
No running production source, geometry, input, optimizer, G3 selection,
holdout label or finite calculation budget is modified.

The existing `BiaxialClampedCurvatureCoordinates` supplies the same controlled
coordinate used by the planned model. In the actual G2 row-cell frame, fix
the open row2 and all fractional atoms while the two prescribed substrate
rows scale as1+epsilon. For this frozen partial derivative, J=dH/d epsilon
has rows0/1 equal to H[:2]/(1+epsilon) and row2 zero. With raw ASE stress,

```text
dE/d epsilon = V * sum(sigma * solve(H, J).T), at P=0.
```

The tensor contraction retains substrate reactions. Projecting them away
would destroy external physical work, though internal relaxation excludes
them. Units are (Angstrom^3)*(eV/Angstrom^3)*(Angstrom/Angstrom)=eV per unit
dimensionless strain. Cell energies divide by4HfO2 then multiply1000 for
meV/HfO2. A1%increment is0.01, not1; no SCF error/uncertainty is estimated.

For comparison, a homogeneous in-plane affine deformation also carries the
open row's tilt: J_affine=H*(I-n*n.T)/(1+epsilon). Its work equals the planar
stress trace. This is a **different frozen internal-coordinate convention**,
not a second relaxed mechanical ensemble or automatically a numerical bug.
The two derivatives differ by raw work on that extra allowed row2 movement.
They agree on a truly stationary open-row branch; an ordinary NEB residual
pass alone does not establish this condition at each sampled image.

Selected raw-derived partial work, meV/HfO2 per unit strain:

| Sample | Fixed open row | Affinely carried open row |
| --- | ---: | ---: |
| PO+ endpoint, image0 | -228.414524 | -228.414524 |
| Sampled highest image3 | -599.344255 | -594.022797 |
| M endpoint, image8 | -2337.156391 | -2336.670108 |

At image3 the difference is5.321459 in these units; the fixed-row work relative
to PO+ is-370.929732. All nine rows/Jacobians/raw stresses are in`final_analysis.json`.
This is a frozen local measurement, **not d(relaxed barrier)/d epsilon**, an
independent forecast, a certified TS or proof of channel-selective protection.
Finite-strain G2 comparison and stationary-branch/strong-control tests remain
required. Numerical work conjugacy here is checked with EMT; it does not
measure native DFT finite-difference stress/energy error or Pulay uncertainty.

Reproduce from the repository root into a new output:

```text
python -m scripts.analyze_hfo2_clamped_external_work --case benchmarks/hfo2_channels/20261008/clamped_M_terminal_E068_20261011 --audit-sha256 56b92e8b15195d6f619a8402c7eee83a774b2e364e012478faa3fdc5d2da431c --output NEW_REPORT.json
```

The study analyzer admits only the current completed M training audit **before
opening labels**. Later training receipts need explicit registration, not an
automatic arbitrary-receipt or held-out-data read. It also checks the original
substrate registration and refuses unknown boundaries/receipt changes; it does
not read or generate the+.005holdout. The reusable physical-work helper remains
independent of this finite-study admission. Tests cover both frozen charts' energy derivatives, tilted/rotated
cells,0/+1%anchors, reaction-only work, nonstationary chart differences and
actual nine-frame replay. Validation distinguishes local and independent
archive tests, source/unit heuristics, native material data and manuscript
reading-layout checks. No final G2/G3 or scientific completion is implied.
