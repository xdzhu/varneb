# Manuscript terminology and claim boundaries

Updated 2026-10-15. Scope: HfO2 manuscript and its source-to-table checks.

| Term | Meaning used here | Does not establish |
|---|---|---|
| Sampled peak | Highest evaluated image of a declared ordinary NEB band | Continuous maximum, stationary TS or barrier error |
| Ordinary-converged band | Meets the declared 0.10 eV/Angstrom NEB residual | Zero physical gradient at each image |
| Internal local curvature | Finite differences in 36 translation-free joint atomic-cell length coordinates at fixed substrate | Mass-weighted phonon frequencies or an uncertainty-resolved stationary index |
| Prescribed control | Substrate strain, excluded from internal instability counting | A free unstable mode |
| Coupled local direction | A direction of the joint atomic-cell matrix in the declared chart | An endpoint phonon eigenmode |
| Gradient-converged candidate | Fixed retained coordinates, released-gradient norm at most 0.03 eV/Angstrom | Stable conditional minimum or unique branch |
| Conditional candidate sheet | Continued candidates at prescribed coordinates, with holes retained | Globally minimized PES or a certified smooth interpolation domain |
| Curvature seed | Standard measured-Hessian quasi-Newton initialization | A new optimizer theorem, guaranteed convergence or whole-path speedup |
| Raw negative count | Negative eigenvalues of a finite-step symmetrized matrix | Physical error bound or certified TS index |
| Censored cost | Calls spent before a common attempt cap without reaching the target | Completed time-to-target or a denominator for a speedup ratio |

Notation is not interchangeable: the representation figures use the declared
mass-weighted atomic T chart; local joint-curvature coordinates have length
units with L=5.12968067458423 Angstrom. Reported E160 total energies and energy
changes use the twelve-atom cell. Existing G2 barrier tables use their stated
per-formula-unit normalization. No ASR, atom rematching, new reference selection
or held-out labels are introduced by the manuscript integration.

The local tests ask whether omitted coordinates are equilibrated well enough
for a reduced description. They do not answer the final independent material
prediction question. Full T optical and all four Cmma registrations remain
strong controls; unequal ranks are disclosed rather than called equivalent.
