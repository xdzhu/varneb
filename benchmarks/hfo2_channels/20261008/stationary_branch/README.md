# E031: fixed-control stationary response, not a material prediction

This completes a missing numerical prerequisite of the registered B2--B5
comparison: after stable R release, follow the retained internal minimum or
first-order saddle at prescribed control t, retaining nonstationary offsets.
The actual mechanical internal basis is explicit and excludes the control.
Unrepresented internal forces and clamped reactions are reported separately.
The algebra is standard; it is not a new theorem or a barrier certificate.
See the [API and units](../../../../docs/CONTROLLED_STATIONARY_BRANCH.md).

## Bounded experiment and actual evidence

-94 focused tests pass in1.93s;38 new cases. One initial integration test
  referenced a nonexistent boundary factory; only that test was corrected to
  the existing `clamped_plane_vcneb_boundary`, then the suite was rerun.
- A clean staged archive excluding user changes passes1116 tests,2skips,
  313warnings in127.54s. No remote pytest or installation is claimed.
- The bounded analytic checker uses8groups and24points: minimum/saddle,
  frozen/stable release and original/rotated charts. Against a separately
  assembled fixed-control solve, clean-local errors are <=2.78e-17 in
  displacement,1.39e-17 in energy and1.25e-16 in represented internal
  gradient. An independent finite-energy control derivative differs by
  <=8.00e-13. These are analytic implementation residuals, not DFT errors.
- HF runs the exact same archived source. All source hashes and83nonfloat
  fields match the clean report;124float fields agree within1.51e-13.
  NumPy versions1.23.5/1.26.4 are retained separately. Every run uses0DFT
  and0calculator calls. The real production source archive remains unchanged.

The saddle control curvature need not soften. In the analytic test,
the moving negative direction changes a frozen control curvature-3 to-2.5,
whereas the corresponding minimum changes it to-3.5. This verifies a sign
distinction in stationary continuation, not a hafnia response or a favorable
switching/escape window. Control-scale and rotated-chart tests preserve the
physical energy while applying the explicit derivative chain rule.

## Retained source-format difference, not a waived check

`analytic_clean.json` is the canonical tested-archive report;
`analytic_HF.json` is its strict same-source remote replay. The original
`analytic_local.json` is also preserved. An initial raw-hash comparison of
that working-tree report failed: the existing `relaxed_curvature.py` has
LF bytes in the worktree and CRLF bytes in the Git archive. Their decoded
lines and newline-normalized bytes are exactly identical. Neither legacy
file nor its raw report hash was rewritten. The current new code/test files
are byte-identical to the clean archive. The delivery receipt retains both
legacy hashes and the failed comparison; archive-versus-HF source checks
remain exact, not weakened or normalized away.

## Reproduction and scientific limits

From the repository root, use a **new** output path:

```sh
bash benchmarks/hfo2_channels/20261008/stationary_branch/autoresearch.sh /new/analytic-check.json
python -m pytest -q tests/test_stationary_branch.py tests/test_quadratic_reduction.py tests/test_relaxed_curvature.py tests/test_biaxial_curvature.py
```

The checker is available with `python -m scripts.check_stationary_branch`
on Windows without Bash. It refuses an existing output before computation.
Full code/data hashes, archive identity and actual test counts are in
[the delivery receipt](validation_delivery.json).

Actual training E/g/H, probe coverage, anharmonicity, full admissible stability,
branch continuity, endpoint-compatible barrier construction and independent
label tests remain required. This does not complete B2--B5 selection/freezing,
G2 material paths, G3 conditional surfaces or JCTC evidence. No new physical
parameter, SCF, production job, optimization restart or held-out label is
introduced to fill waiting time. Both existing G1 switching handles continue
under the original finite budget and ordinary0.10 criterion.
