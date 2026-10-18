# One bounded preserving-path connection test

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
are retained. No physical parameter changes or sign/strain search.

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
actual runtime admission. No convergence or new barrier is yet claimed.

Job28949435 started onnode1 and reproduced the cached step-zero residual
0.125611 without new endpoint calculations. Its first fresh interior SCF
completed in85.041s with genuine32MPI, unchanged six input assets and full
finite E/F/stress; see `startup_snapshot.json`. This is one completed
native call, not a complete NEB update or convergence claim. Only E269
preserving and this connection job occupy study allocations.

```text
python benchmarks/hfo2_channels/20261008/preserving_connection_E273_20261019/prepare.py --output FRESH_DIRECTORY
python -m pytest -q tests/test_hfo2_preserving_native_connection.py tests/test_hfo2_preserving_connection_seed.py
```
