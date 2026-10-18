"""Ordered atom-level E273 comparison; geometry only, no DFT or classification fit."""
from __future__ import annotations
import argparse
import csv
import importlib.util
import json
from pathlib import Path
import tempfile

import numpy as np
from ase.io import read
from scripts.audit_hfo2_static_replica import sha256
from scripts.analyze_hfo2_clamped_reference import original_chart_images
from scripts.plot_hfo2_preserving_connection import restore, native_images

REPO = Path(__file__).resolve().parents[1]
BASE = REPO/'benchmarks/hfo2_channels/20261008'
CASE = BASE/'preserving_connection_E273_20261019'


def band_geometry(images, reference, orientation):
    """Use the existing proper rotation/row order, never atom permutation or MIC."""
    images = original_chart_images(images, orientation, reference)
    frac = np.array([a.get_scaled_positions(wrap=False) for a in images])
    if frac.shape != (9, 12, 3): raise ValueError('registered nine-image Hf4O8 band required')
    displacement = (frac-frac[0])@reference.cell.array
    relative = displacement-displacement[:, :4].mean(axis=1, keepdims=True)
    # Exact registered atom identities in the original T chart. These are
    # geometrical layer bands, not a space-group or literature-family assignment.
    lower_pair, upper_pair = [0, 3], [1, 2]
    z = frac[:, :, 2]
    margins = []
    for at in range(4, 12):
        lo, hi, lift = (lower_pair, upper_pair, 0) if z[0, at] < .5 else (upper_pair, lower_pair, 1)
        margins.append(np.minimum(z[:, at]-z[:, lo].max(axis=1),
                                 z[:, hi].min(axis=1)+lift-z[:, at]))
    margins = np.array(margins).T
    if not np.isfinite(margins).all(): raise ValueError('finite ordered layer geometry required')
    return dict(fractional_positions=frac.tolist(), Hf_relative_displacement_A=relative.tolist(),
        minimum_layer_band_margin_fractional_z=float(margins.min()),
        layer_band_margins_fractional_z=margins.tolist(),
        no_Hf_layer_band_crossing_at_samples=bool((margins > 0).all()),
        max_Hf_layer_corrugation_fractional_z=float(max(np.ptp(z[:, lower_pair], axis=1).max(),
                                                       np.ptp(z[:, upper_pair], axis=1).max())),
        center_O_mean_relative_y_A=float(relative[4, 4:, 1].mean()),
        center_O_max_abs_relative_y_A=float(np.abs(relative[4, 4:, 1]).max()),
        center_O_relative_y_nonuniform_norm_A=float(np.linalg.norm(
            relative[4, 4:, 1]-relative[4, 4:, 1].mean())))


def analyze():
    path = BASE/'preserving_pattern_branches_E272_20261018/analyze_prefix.py'
    spec = importlib.util.spec_from_file_location('E273_atom_motion_original', path)
    original = importlib.util.module_from_spec(spec); spec.loader.exec_module(original)
    old, _ = original.original_images()
    reference = read(BASE/'reference_variants/T.vasp')
    orientation_path = BASE/'clamped_endpoint_seeds/clamped_seed_manifest.json'
    orientation = json.loads(orientation_path.read_bytes())['orientation']
    with tempfile.TemporaryDirectory(prefix='E273-motion-') as td:
        folder, outer = restore(9, Path(td), terminal=True)
        final, _, receipt = native_images(folder, 9)
        bands = {name: band_geometry(images, reference, orientation)
                 for name, images in [('original', old), ('lower', final)]}
        native_source = dict(archive_sha256=outer['archive_sha256'], source_job=receipt['job_id'], step=9)
    a, b = (np.array(bands[name]['fractional_positions']) for name in ('original', 'lower'))
    delta = (b[4]-a[4])@reference.cell.array
    delta -= delta.mean(axis=0)
    hf_relative = delta-delta[:4].mean(axis=0)
    return dict(status='ordered_geometry_comparison_not_mode_energy_or_TS_certificate',
        source=native_source, original_observation_sha256=original.OBSERVATION_SHA,
        orientation_sha256=sha256(orientation_path), reference_sha256=sha256(BASE/'reference_variants/T.vasp'),
        analysis_source_sha256=sha256(Path(__file__)), atom_order=['Hf']*4+['O']*8,
        coordinate_frame='existing original-T chart; fixed T metric; band motions measured from each common initial endpoint',
        bands=bands, center_lower_minus_original_translation_free_A=delta.tolist(),
        center_lower_minus_original_Hf_relative_A=hf_relative.tolist(),
        center_difference_norm_A=float(np.linalg.norm(delta)),
        oxygen_fraction_squared_difference_norm=float(np.sum(delta[4:]**2)/np.sum(delta**2)),
        new_DFT_calls=0, literature_CN_ST_family_assigned=False,
        electronic_polarization_computed=False, phonon_eigenmode_assigned=False,
        limitations=['Layer-band margins are ordered fractional-z geometry, not physical distances.',
            'No crossing at the nine samples; positivity also holds for their piecewise-linear fractional interpolation.',
            'No statement about an unknown continuous relaxed path, TS index, irrep or electronic polarization.',
            'Nonzero individual lateral motions can have a zero species-mean transverse proxy.',
            'The old/new centers do not hold all observables fixed; matched causal control is the separate E272 F/A/J test.'])


def write_report(output):
    report = analyze(); output = Path(output); output.mkdir(parents=True, exist_ok=False)
    with (output/'analysis.json').open('x', encoding='utf-8', newline='\n') as h:
        json.dump(report, h, indent=2, allow_nan=False); h.write('\n')
    with (output/'center_motion.csv').open('x', encoding='utf-8', newline='') as h:
        writer = csv.writer(h); writer.writerow(['atom_index_zero_based', 'species', 'dx_A', 'dy_A', 'dz_A',
                                                 'Hf_relative_dx_A', 'Hf_relative_dy_A', 'Hf_relative_dz_A'])
        for i, symbol in enumerate(report['atom_order']):
            writer.writerow([i, symbol]+report['center_lower_minus_original_translation_free_A'][i]
                            +report['center_lower_minus_original_Hf_relative_A'][i])
    return dict(output=str(output), new_DFT_calls=0,
                center_difference_norm_A=report['center_difference_norm_A'])


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    print(json.dumps(write_report(p.parse_args().output)))
