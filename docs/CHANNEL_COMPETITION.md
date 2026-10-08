# Same-initial-state channel competition

`vcneb.ChannelPath`, `summarize_competing_paths` and
`compare_channel_selectivity` consume numeric, already evaluated paths. They
do not call a calculator, optimize a path, alter atom identities, or certify
a stationary saddle. The implementation is an auditable analysis contract,
not a new kinetic theory or a novelty claim by itself.

## One physical reference and one declared set of candidates

Every path must have the same ordered periodic initial structure, raw initial
energy, physical-input hashes, mechanical family/parameters, and pressure.
An integer lattice shift common to a complete chain is acceptable; atom
permutations and independent image-wise relifting are not. A reversed view
of T→PO supplies PO→T without a new DFT calculation. It retains the original
source optimization metric and residual rather than recomputing a nominally
different NEB metric from the reversed endpoint.

For raw cell energies E_i and volumes V_i, the diagnostic uses
H_i=E_i+P V_i, with P in eV/A³. It reports both source directions and the
direction from the common initial state, in eV per formula unit. Including
endpoints in the discrete maximum ensures nonnegative directional barriers.
For each source, B_forward−B_reverse=H_final−H_initial; a reversed display
changes the sign of its reaction energy, not the raw data.

Pass `required_channels={name: "switching" or "decay", ...}` explicitly.
The minimum for a role is absent when any required path in that role is
missing. Unknown/duplicate channels are rejected rather than silently
changing the declared set. This guards against overlooking PO→T while only
showing a favourable comparison against PO→M. It does **not** prove that all
physically possible channels have been enumerated.

## Provisional maxima are not measured-error barriers

The module always exposes valid raw discrete profiles, including unconverged
ones. `ready_for_bounded_discrete_comparison` remains false until every
required path passes the ordinary residual threshold, has a sampling-audit
SHA256, and supplies a nonnegative measured barrier-error bound together
with its audit SHA256. The caller must actually audit those linked data:
the numeric API records provenance but does not inspect an arbitrary
external audit document. A 0.10-eV/A residual is never converted into an
energy uncertainty. A missing error is never treated as zero.

Error bounds apply to **directional barriers**, in eV per simulation cell,
not separately fitted endpoint errors or Gaussian standard deviations. The
adapter for the current HfO₂ observations supplies no such bounds because
none have yet been measured. Accordingly, its selectivity verdict is gated.

## Separate relative selectivity from absolute barrier changes

For the declared candidates define B_s=min(B_switch), B_d=min(B_decay),
and S=B_d−B_s. If each barrier lies in [L_j,U_j], its minimum lies in
[min(L_j),min(U_j)]; the identity of the lowest channel may change inside
the uncertainty range. Difference intervals use [L_left−U_right,
U_left−L_right]. This conservative arithmetic does not assume independent
errors and may overestimate uncertainty when sources are correlated.

`compare_channel_selectivity` requires the same physical inputs, declared
channel set, formula normalization, pressure, and mechanical family. It
allows parameter values such as strain to change, but not their definitions.
Released zero-pressure and zero-strain clamped results are separate boundary
families, not a smooth strain derivative.

A positive lower bound for ΔS supports increased **relative selectivity**.
It does not establish increased resistance to decay: both B_s and B_d can
fall while S rises. The stronger diagnostic requires an upper bound below
zero for ΔB_s and a lower bound for ΔB_d no smaller than −δ, where the
noninferiority margin δ must be declared before inspecting results.
Error-bar overlap alone does not prove noninferiority. The API requires δ
explicitly; it neither chooses a physical significance threshold nor proves
that the caller preregistered it. Setting δ=0 is the strict version.

No small-cell, zero-temperature barrier diagnostic establishes switching
rates, retention time, coercive field, finite-temperature stability, or the
absence of unexamined leakage paths.

Units remain explicit in the numeric ASE API: raw E in eV/cell, pressure
in eV/A³, cell vectors in A, barrier errors in eV/cell, and outputs in eV/fu.
Regression independently checks pressure-volume work against SI joules and
the library electron charge. The standard-library unit/uncertainty static
auditor reports no findings in the two new analysis modules; this heuristic
check is not a complete dimensional proof. Pint is absent from the current
production-compatible Python environment; no dependency or physics constant
is changed for this audit. Deterministic numerical bounds are not Gaussian
standard uncertainties or 95% confidence intervals.

## Reproducible HfO₂ example

The first four-channel dataset reuses 37 complete SCF image evaluations:
T→PO step10 viewed in reverse, PO→M step10, and both flip candidates step1.
All share PO+, P=0/E=0, ABACUS100Ry/full10auDZP. Their unequal iteration
stages and unmeasured errors prevent a final channel-ranking claim.

```bash
python -m scripts.analyze_hfo2_channel_network \
  --specification benchmarks/hfo2_channels/20261008/channel_network/network_specification.json \
  --output /a/new/path/channel_analysis.json
python -m pytest -q tests/test_channel_competition.py tests/test_hfo2_channel_network.py
```

The output must be fresh. The script validates numeric E/F/stress, POSCAR
and trajectory hashes, physical contracts, source optimizer rows and force
replay before comparison. Continuous T-reference projections remain
descriptive, not local TS phonons or modal energy partitions. See
`benchmarks/hfo2_channels/20261008/channel_network/README.md` for observations
and boundaries; source raw SCFs and all live calculations remain untouched.
