# Physical force diagnostics under mechanical constraints

`VCNEB.path_diagnostics()` and `VCNEB.saddle_diagnostics()` distinguish the
ordinary NEB residual from the physical force. The residual removes a physical
path-parallel force and adds springs, so a residual pass alone is not a TS
certificate. The default ordinary threshold remains0.10eV/Angstrom.

Under a prescribed mechanical or mode constraint, distinguish a physical
gradient in the allowed space from a support reaction. With the active-force
mask already applied, the implemented orthonormal subspace projector P gives
`F_allowed = P @ F_active`, and `F_reaction = F_active - F_allowed`.
Masks and the registered generalized-coordinate metric remain part of this
definition; a force already removed by a coordinate mask is not reconstructed.
This is the actual optimizer constraint, not a basis selected merely for
plotting or an omitted but physically movable response direction.

Both public summaries include these explicitly named fields:

| Field suffix, in eV/Angstrom | Meaning |
|---|---|
| `allowed_true_generalized_force_norm_eV_per_A` | Euclidean norm of the physical allowed force |
| `allowed_true_generalized_force_max_vector_eV_per_A` | Maximum three-component allowed force norm, under the original metric |
| `constraint_reaction_force_norm_eV_per_A` | Euclidean norm of the projector-excluded active force |
| `constraint_reaction_force_max_vector_eV_per_A` | Maximum three-component norm of that excluded force |

The saddle summary additionally includes `allowed_true_tangential_force_eV_per_A`
and `allowed_true_perpendicular_force_eV_per_A`. Their meanings match the
already projected `true_tangential_force_eV_per_A` and
`true_perpendicular_force_eV_per_A` in the per-image **path** summary.

For compatibility, legacy `true_*` norms in the **saddle** summary retain the
active raw force before the mode/subspace projector. They can contain reactions
and must not be used as an allowed-space convergence test. The raw Cartesian
stress and existing path fields are still separately available. With no
subspace projector, new allowed and old active-raw norms agree; the
projector-excluded active reaction is zero.

The same physical gradients may have different max-vector and Euclidean norms;
neither is a kbar traction, a barrier-energy error estimate or a Hessian index.
Endpoint traction screens retain their own declared units and thresholds.
Excluding a support reaction from **internal stationarity** does not set its
physical value to zero. Responses to a change of the prescribed substrate
must still use the raw stress and the actual cell Jacobian/volume; dropping
the reaction there would remove part of the external strain-work derivative.
Neither removing a prescribed reaction nor satisfying the ordinary NEB target
proves a continuous index-one saddle or permits deleting an unstable mode.

The [native HfO2 cache example](../benchmarks/hfo2_channels/20261008/clamped_stationarity_E073_20261011/README.md)
reproduces these distinctions without new DFT and with unchanged ordinary force.
