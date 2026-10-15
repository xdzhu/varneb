# Atomic versus joint release: first measured matched treatment comparison

2026-10-15. This consumes an actual native audit, not synthetic atomic
labels. It does not alter the calculator or fit the original predictions.

## The physical question and current answer

在固定两个路径坐标后，原先的大幅降能主要来自原子重新排布，还是额外的
晶格放松？目前翻转路径中心点的实测答案是：只放松原子降低48.348691
meV/胞，再允许既定晶格自由度与原子联合放松，额外降低1.136239 meV/胞。
因此，不能把该点的大幅降能主要归因于晶格放松。这是指定处理之间的差，
不是唯一的模式能量分解，更不能推广为晶格对整条路径不重要。

All energies below are **meV per twelve-atom Hf4O8 cell**. Divide by four
for meV/HfO2. F freezes omitted coordinates; A releases the31atomic fixed-Q
directions at that Q seed's fixed cell; J releases both atomic and permitted
cell coordinates. Different Q seeds need not have the same cell.

| Existing raw center, Q=(0,0) | A-F | J-A | J-F |
|---|---:|---:|---:|
| T escape | -0.271398 | -0.380744 | -0.652142 |
| Preserving flip | -48.348691 | -1.136239 | -49.484930 |

These are local candidates, not certified global nested minima, TSs,
stable conditional branches, full-path barriers or independent predictions.
In particular the atomic-gradient passes do **not** imply joint stationarity:
the joint released norms at the two atomic origins are0.110430/0.194104eV/A,
both above the separate0.03released-gradient target.

## Complete T axial cohort, incomplete preserving cohort

The fresh E176 native prefix reparsed20new SCFs and exactly replayed all
seven completed node optimizations. It covers T0--4 and preserving0--1.
The original job was RUNNING; preserving2 was absent when the prefix was
read. There are no missing T axial labels, but only **one of four** measured
non-origin preserving atomic labels. No partial preserving cohort is ranked.

Both original probe-step models and all20model-point rows are retained.
For T's same four axial points, origin-aligned maximum model discrepancies
are:

| Original Hessian step/A | F | A | J |
|---|---:|---:|---:|
| 0.01 | 0.416264 | 0.501879 | 1.327177 |
| 0.02 | 0.403160 | 0.515704 | 1.346071 |

The differences between these maxima are **not** a physical decomposition:
the corresponding error identities are evaluated pointwise in analysis.json.
The larger T joint shape discrepancy warrants checking the additional
release response, not declaring a microscopic mechanism or a physical error
bound from these small differences. Preserving's large origin release is
mostly present in A, while its remaining three atomic axial labels and
released stability are not yet audited here. No preferred step is selected.

## Provenance, reproduction and costs

- atomic_native_prefix.json:673714cb35951ee60b856a1bfb7775aba91ab6ee8a3c2cd60bbe7d8c50dd2e96.
- analysis.json:16b9dde68203a9888a642959c2a6ce979cd448f9ceb926a0cec911d06d61d5b5.
- The original E176 pre-label20-row forecast and E178 postprocessor are
  unchanged. The native reader uses original E161 runtime and original
  E175 worker, not this checkout as production source.

From the repository root, with numerical threads1:

```text
python benchmarks/hfo2_channels/20261008/G3_treatment_comparison_E178_20261015/analyze.py --atomic-audit benchmarks/hfo2_channels/20261008/G3_atomic_treatment_E192_20261015/atomic_native_prefix.json --atomic-audit-sha256 673714cb35951ee60b856a1bfb7775aba91ab6ee8a3c2cd60bbe7d8c50dd2e96 --output NEW_COMPARISON.json
```

The supplied prefix covers20new atomic SCFs:9atT,11atpreserving, with
22.331631completed-call proxy core-hours. In-flight or later records are
outside that audit, not free. Eight old frozen statics/8.136814proxy
core-hours, original296curvature SCFs, all joint paths/preparation and the
currently running released-curvature batch are additional. This is not a
total campaign cost or speedup comparison.

No new DFT was initiated, scheduler state changed, production source edited,
parameter retuned, independent label consumed, or physical+.005holdout read.
The original five user edits and main manuscript are preserved. The next
physical action remains the full B3 and released-curvature audits, followed
by the registered admission/forecast checks, not an enlarged calculation.

Related verification:75local and75independent checks across four modules,
including10new checks, zero failures/errors/skips. This is not the entire
repository suite. The independent CLI reproduces analysis.json bytewise.
Tested source tree3fb3893bc83a7503cf19787bad4b501e1eab24a8/archive
809400639dc71c3d0c62591a1bd5d91ec5ffedfc7dc22ca3396cc85cbe8f00db.
All scratch, test temporaries and independent exports are on C:.
