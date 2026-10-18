# One bounded preserving-path connection test

**Completed28949435, 2026-10-19 at02:27:41 CST.** Nine ordinary updates
reached0.0940120033eV/A, using63fresh interior SCFs and no new endpoint
calls. The independently replayed sampled peak is83.1662809401meV/HfO2,
13.6619020627 below the original96.8281830028. The central transverse
structural proxy remains0.1210744444Angstrom. This lower connected candidate
changes the interpretation of the original path set: its1.250meV/HfO2
separation from the original T escape is not a resolved ordering at the
current finite sampling and ordinary target. Neither a stationary TS nor
a global MEP is certified, and old forecast labels/scores are unchanged.

The [final mechanism figure](figures/terminal_step_0009/README.md) links
energy, transverse/z structure and all interior residuals. Its source
reader independently reparses all63native calls, not only the last frame.
The519-file complete archive `native_terminal.tar.gz` has SHA256
`54bda58fa6ea2e2bf1bbc6ae1ef193379ebb90aa5a386f699e6546b0e6056d10`.
New costs are55.718818native-call proxy/55.857778allocation core-hours;
earlier preparation/release/path costs remain additional. Full pseudo/orbital
bytes were rechecked on HF, not redistributed or rehashed locally.

The records below retain the original admission and progress stages.

**Submitted and released as28949435 after native preflight.** The original
E272 job completed0:0, freeing the second study slot. All nine exact
native caches and structure roundtrips passed with zero fresh DFT; the
replayed initial NEB residual is0.1256114037eV/A. Held one-node32CPU/6h
resources were verified before release. Runtime source3f72b35f was pushed
before deployment; the ten-file immutable bundle pin is
`a452d5b125d6d95c729ea60607524ce5d35370d5aab57ac6f2b1b7869d61991a`.
See deployment/native-preflight/held/accepted/released receipts. No
original E269/E272 runtime source was overwritten or third job submitted.

The [protocol](protocol.md) asks whether the lower, matched-observable
zero-strain J branch can form a connected ordinary-NEB path. The current
cached polygon has85.499889meV/HfO2 sampled maximum but0.125611eV/A NEB
residual, so neither a revised barrier nor a TS is claimed. All seven
interiors evolve freely; the local intervention constraints do not become
new NEB constraints. The original two endpoints and atom/periodic gauge
are retained. This was the pre-relaxation seed, not the final result above.
No physical parameter changes or sign/strain search.

The first complete free-chain update at00:55:24 has been exported and
independently replayed: sampled maximum85.269291meV/HfO2, residual0.123740,
still unconverged. The [four-panel mechanism preview](figures/step_0001/README.md)
shows original/seed/step1 energies, transverse/z structural proxies and
per-image residuals. Signed component checks locate the remaining force
primarily in the perpendicular physical term, not spring spacing. No k
change, new prediction score or acceleration claim follows. The76-file
selected snapshot stays compressed in `snapshot_step_0001.tar.gz`; it
contains only this evaluated frame, not all past/future runtime calls.

The inert `prepared/` seed includes all nine exact native E/F/stress cache
references, the single replaced J geometry, the full connection replay
and byte hashes. Standard public resume-cache behavior invalidates each
interior's results after it moves. No new SCF or scheduler mutation occurs
in preparation. A real run is bounded by24updates/168interiorSCFs/6h,
with the original0.10eV/A target and no CI.

`control.py`, `run.slurm` and `package.py` implement the finite entry point
without submitting anything. Native preflight refuses a live/timeout parent,
checks all nine exact E/F/stress caches with fresh DFT disabled and replays
the same initial NEB residual. Execution requires the exact separately
accepted one-node32CPU/6h identity, one thread/rank and at most two study
allocations. The deployed runtime remains the original immutable E268
source, not a mutable local checkout. Preparation metadata retains its
original non-admission status; the later native/held receipts establish
actual runtime admission. This describes admission, not the later terminal result.

Job28949435 started onnode1 and reproduced the cached step-zero residual
0.125611 without new endpoint calculations. Its first fresh interior SCF
completed in85.041s with genuine32MPI, unchanged six input assets and full
finite E/F/stress; see `startup_snapshot.json`. This is one completed
native call, not a complete NEB update or convergence claim. Only E269
preserving and this connection job occupy study allocations.

`collect_terminal.py` is a separate read-only postprocessor, not a runtime
replacement. It requires the original COMPLETED0:0 handle, distinguishes
native0.10 convergence from the24-update step cap, and retains every fresh
SCF and its cost. The nine initial caches mean exactly7*updates fresh
interior evaluations and no fresh endpoint calls. It neither relaunches
a bounded failure nor promotes a live/timeout/exit-zero-only result.
It was executed only after the original handle became COMPLETED0:0.

```text
python benchmarks/hfo2_channels/20261008/preserving_connection_E273_20261019/prepare.py --output FRESH_DIRECTORY
python -m pytest -q tests/test_hfo2_preserving_native_connection.py tests/test_hfo2_preserving_connection_seed.py
python -m pytest -q tests/test_hfo2_preserving_connection_terminal.py
```
