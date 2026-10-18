# One bounded preserving-path connection test

**Prepared only, not submitted.** E272 must finish and one of the two study
allocation slots must be free before a separately frozen runtime admission.
The original E269/E272 running sources are unchanged.

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
source, not a mutable local checkout. These guards are not a claim that
preflight or a production run has already succeeded.

```text
python benchmarks/hfo2_channels/20261008/preserving_connection_E273_20261019/prepare.py --output FRESH_DIRECTORY
python -m pytest -q tests/test_hfo2_preserving_native_connection.py tests/test_hfo2_preserving_connection_seed.py
```
