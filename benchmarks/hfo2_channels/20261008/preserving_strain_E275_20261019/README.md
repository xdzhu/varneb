# Strain response of the lower connected preserving candidate

Preparation only; no job or result is claimed by this document.

The finite [protocol](protocol.md) transfers the actual E273 step9 lower
zero-strain candidate to the already studied +0.5% and +1% boundaries,
reusing their original fixed endpoint caches. The two complete free bands
ask whether transverse atomic accommodation changes the strain-dependent
competition with T escape. This is not another forecast competition.

Two one-node32CPU allocations at most, ordinary0.10/noCI,60updates/18h
per band,427fresh interior calls maximum per band. All original electronic
settings and data labels remain unchanged. No automatic extra study follows.

```text
python benchmarks/hfo2_channels/20261008/preserving_strain_E275_20261019/prepare.py --output FRESH_DIRECTORY
python -m pytest -q tests/test_hfo2_preserving_strain.py
```
