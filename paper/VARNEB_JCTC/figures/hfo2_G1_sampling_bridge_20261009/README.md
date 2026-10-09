# Six actual sampling checks on the preserving G1 reconstruction

Figure conclusion: real same-contract static points exceed the original
nine-image sampled maximum on the specified fractional/cell reconstruction.
An ordinary residual pass does not determine this between-image maximum.
This is not an independently relaxed continuous MEP, stationary TS, error
bound or prospective prediction benchmark.

Panel(a) places nine reused source samples and six new SCFs on the original
normalized generalized arc. Panel(b) shows each tested segment at fractions
0,0.25,0.50,0.75,1. Straight connections are sample guides; the dotted line
marks the old sampled maximum32.806023meV/HfO2. New highest observed
sample38.508448 exceeds it by5.702425meV/HfO2. Both sides are separately
evaluated, not imposed as mirrors.

The [material case](../../../../benchmarks/hfo2_channels/20261008/preserving_sampling_bridge_20261009/README.md)
contains the fixed six-point recipe, actual job28374431 provenance, all
raw E/F/stress and input/log/geometry audits. The old nine-image step69
is unchanged. All fifteen records replay before plotting; the inserted
cached band has0.099372458eV/A ordinary residual with zero optimizer steps
or further SCFs. There are no independent replicates or statistical bars.

Exports: editable-text SVG and PDF; 600-dpi LZW TIFF; PNG preview;
`source_data.csv` with all fifteen rows; `qa.json` with source/hash metadata.
Dimensions183x99mm, aligned shared-energy axes, plain larger(a)/(b), no
panel titles, all spines visible and translucent framed white legends
outside the data area. Actual rendered-pixel review is recorded in the
delivery receipt; the reproducible generated QA retains its pending flag.

Reproduce from repository root:

```powershell
python -m scripts.plot_hfo2_sampling_bridge --observation benchmarks/hfo2_channels/20261008/switching_converged_update_20261009_1255/preserving_step69 --completed benchmarks/hfo2_channels/20261008/preserving_sampling_bridge_20261009/completed_HF --output NEW_FIGURE_BUNDLE
```

The output must not exist. The script cannot execute DFT or start another
band. It validates real raw logs, full energy/forces/stress, physical
input hash records, source geometry/writer equivalence and actual energy
conversion before drawing. Large repeated orbital/pseudopotential payloads
are not duplicated in the repository; their original bytes are retained
and audited on HF. No retuning of100Ry/full10auDZP or other INPUT/KPT
settings occurred.
