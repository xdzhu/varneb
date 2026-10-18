# Complete original-target analysis

`analysis.json` replays ten actual bands and90retained native image rows;
`outcomes.csv` records three same-strain comparisons; `profiles.csv`
retains every plotted original profile, including the upper-initialized
midpoint. All42fresh target calls are independently checked against
their native logs, E/F/stress results, geometry and original inputs.
Full pseudopotential/orbital byte checks occurred onHF, not locally.

This final analysis supersedes the *status*, not the bytes, of the earlier
`completed_analysis_0050` snapshot. Original training/independent labels
and forecasts are not replaced. No new independent prediction, global
MEP/TS certificate or energy error bound is inferred. Six preserving
peak geometry checks retain all three symmetry tolerances, software
versions and warnings; no coordinates were standardized or reordered.

New target cost is42SCFs,35.798164call-proxy/36.008889allocation
core-hours. Earlier source/release/curvature/path costs are additional.
The archive is not charged again when read or plotted.

```text
python benchmarks/hfo2_channels/20261008/preserving_strain_E275_20261019/analyze_terminal.py --output FRESH_DIRECTORY
```
