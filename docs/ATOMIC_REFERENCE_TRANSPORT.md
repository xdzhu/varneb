# Atomic reference spans in a local response chart

`transport_atomic_reference_span` is a backend-independent tangent/unit/gauge
adapter, not an eigensolver, Hessian substitute or prediction algorithm.
Use it when a registered atomic reference is mass-weighted in one row-cell
chart but the measured local response uses unweighted lengths in another.
No calculator is invoked and no atoms, forces or source vectors are modified.

## Declare the actual correspondence

Inputs are ordered source ASE atoms, a previously declared **real** mass-
weighted atomic span V, its ordered masses, a supported local response chart
and the registered unimodular integer cell-basis matrix P. Complex vectors
must first pass the separate real-subspace policy; imaginary parts are not
discarded. The source and target may differ physically, but their registration
is the caller's responsibility. No atom assignment, fitted frame, MIC or
integer relift is guessed by this API.

At the target chart origin, with ASE row cells, the raw atomic tangent is

`A = (V/sqrt(m)) @ inv(H_source) @ inv(P) @ H_target`.

The multiplication acts on each atom's three Cartesian components. P is a
cell-basis change, not interchangeable with a Cartesian rotation. The target
cell already encodes its Cartesian frame. All cell and prescribed-control
components of A are zero: this does not freeze the cell when a later solver
releases its separately declared allowed complement.

The chart supplies three translation directions, which are normalized for
the analysis projector T. Let `A_internal=(I-T@T.T)@A`. The returned B has
orthonormal columns and C is the source-amplitude map:

`A = B@C + A_translation`, `q_local=C@Q_source`.

Q_source has sqrt(amu)Angstrom units; q_local has Angstrom units and C has
1/sqrt(amu) units. SVD axes are mixed geometric directions, not named phonons
or new frequencies. Any loss of declared rank is refused. The algebraic rank
tolerance is not a curvature-resolution floor or a physical energy error.

## Keep the excluded work, not an ASR repair

For raw local covector g and source-amplitude increment dQ,

`g.T@A@dQ = g.T@B@C@dQ + g.T@A_translation@dQ`.

No raw net force is zeroed. Restricting an analysis gauge is not a claim that
finite-basis DFT is perfectly translation invariant; the discarded term is
reported. Frozen source arrays stay unchanged. Orthonormalization changes
coordinate amplitudes: never copy source frequencies, LDA force constants or
an old Hessian into these length coordinates without a consistent physical
measurement and covariant transformation.

Pass B to the existing measured-column response controls only after proving
its coverage in the actual measured span. Full cell/control Hessian columns,
stable release, appropriate internal index, branch connections and domain
coverage remain separate requirements. Complete valid harmonic release is
still invariant to a reference rename; this adapter creates no accuracy gain.

## Material example and data scope

The [E154 HfO2 example](../benchmarks/hfo2_channels/20261008/reference_transport_E154_20261014/README.md)
uses the two already frozen E150 training centers and original T/Cmma sources.
It retains ranks3/4/33 and all four separate Cmma registrations. The full
optical T bank is available from an existing matrix, not a newly submitted
phonon calculation. Source preparation and shared curvature costs remain
chargeable; using a cached input is not a zero-cost complete model.
Only scalar checks/source hashes are committed, never private author arrays.

This step neither selects conditional axes nor opens a held-out geometry,
enlarges the finite material matrix or changes ABACUS100Ry/full10auDZP.
