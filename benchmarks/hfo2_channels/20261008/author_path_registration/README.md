# PO-to-M published seed registration

Source: Ma and Liu, [PRL130096801 (2023)](https://doi.org/10.1103/PhysRevLett.130.096801),
[public atomic structures at commit2804a091](https://github.com/sliutheorygroup/supplementary-material/tree/2804a091ad7d36d3f6d068b9e58da5139ce9f8e8/L23),
`figure2/figure2a/transition-poscars/pca21-M`.

`published_seed.traj` stores the twenty parsed structure records already used
by the existing VARNEB case. `registered_author_path.traj` applies the single
declared reflection/translation/permutation in `registration.json` to all
twenty records. This is a coordinate-gauge registration, not recalculated DFT.
The target PO structure is the cached same-contract ABACUS PO+ endpoint.
The small residual functional/relaxation difference at the initial endpoint
will be removed by endpoint correction before a new production chain.

External atomic-coordinate source attribution is retained; the software
license does not assert relicensing of external data. No VASP POTCAR,
wavefunction, private electronic input or third-party manuscript is included.
