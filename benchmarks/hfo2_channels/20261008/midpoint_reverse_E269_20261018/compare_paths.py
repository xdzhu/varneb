"""Compare completed midpoint paths without replacing either original label.

Existing Python/matplotlib paper workflow. Evidence contract: same target
endpoints and boundary; two initial polygons; energy, oxygen rearrangement
and open-axis response shown together. No smoothing or barrier error bars.
"""
import argparse
import csv
from io import StringIO
import json
from pathlib import Path
import sys

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[4]))

import numpy as np
from ase.io import read

from examples.hfo2_fixed_input_factory import CONTRACT, same_ordered_geometry
from scripts.analyze_hfo2_physical_holdout import analyze, BASE, RESULTS
from scripts.audit_hfo2_static_replica import audited_results, sha256
from scripts.plot_hfo2_mechanism import (
    finish, original_chart_images, path_coordinates, rotated_t_triplet, style,
)

CASE = BASE/'midpoint_reverse_E269_20261018'
PIN = '381ae9787cca91b1435b121dafd74ee072d1aa3bd22335f05a7e7ba42bfc7ae1'
JOBS = {'preserving': '28946180', 'T': '28946181'}


def paired_summary(rows):
    if len(rows) != 2 or {r['channel'] for r in rows} != {'preserving', 'T'}:
        raise ValueError('both fixed channels required')
    by = {r['channel']: r for r in rows}
    required = [r[k] for r in rows for k in ('from_zero_peak_meV_fu', 'from_upper_peak_meV_fu')]
    if not np.isfinite(required).all() or min(required) < 0:
        raise ValueError('finite nonnegative endpoint-inclusive peaks required')
    old_gap = by['T']['from_zero_peak_meV_fu']-by['preserving']['from_zero_peak_meV_fu']
    new_gap = by['T']['from_upper_peak_meV_fu']-by['preserving']['from_upper_peak_meV_fu']
    return dict(from_zero_T_minus_preserving_meV_fu=old_gap,
        from_upper_T_minus_preserving_meV_fu=new_gap,
        seed_gap_difference_meV_fu=new_gap-old_gap,
        per_channel_upper_minus_zero_meV_fu={c: r['from_upper_peak_meV_fu']-r['from_zero_peak_meV_fu'] for c, r in by.items()},
        physical_barrier_error_bound=None, original_prediction_scores_replaced=False,
        complete_four_channel_ordering_claimed=False)


def compare(repository, exported, *, initial=False, channels=('preserving', 'T')):
    repository, exported = Path(repository), Path(exported)
    channels = tuple(channels)
    if not channels or len(set(channels)) != len(channels) or not set(channels) <= set(JOBS):
        raise ValueError('explicit unique original channels required')
    original = analyze(repository)
    t = read(repository/BASE/'reference_variants/T.vasp', format='vasp')
    parent, basis, _ = rotated_t_triplet(t)
    orientation = json.loads((repository/BASE/'clamped_endpoint_seeds/clamped_seed_manifest.json').read_bytes())['orientation']
    bands, outcomes = [], []
    for channel in channels:
        old_band = next(b for b in original['bands'] if b['channel'] == channel and b['strain'] == .005)
        old_audit = json.loads((repository/RESULTS/channel/'audit.json').read_bytes())
        old_obs = json.loads((repository/RESULTS/channel/f"step_{old_audit['step']:04d}"/'observation.json').read_bytes())
        folder = exported/channel
        if not initial and not (folder/'audit.json').exists() and (folder/'initial_receipt.json').exists():
            raise ValueError('initial snapshots are not completed paths; request --initial explicitly')
        audit = json.loads((folder/('initial_receipt.json' if initial else 'audit.json')).read_bytes())
        source = folder/f"step_{audit['step']:04d}"/'observation.json'
        if (audit['job_id'] != JOBS[channel] or audit['bundle_sha256'] != PIN
                or sha256(source) != audit['observation_sha256']):
            raise ValueError('original pinned E269 native observation required; never relabel E268')
        if initial:
            if (audit['status'] != 'complete_initial_evaluation_not_terminal_result'
                    or audit['terminal_path_claimed'] or audit['step'] != 0):
                raise ValueError('only explicitly nonterminal initial snapshots are accepted here')
        elif (audit['status'] != 'ordinary_reverse_midpoint_audited' or not audit['ordinary_converged']
                or not audit['all_native_calls_reparsed'] or not audit['genuine_MPI32_checked']):
            raise ValueError('completed ordinary native results required; initial snapshots are not final paths')
        observation = json.loads(source.read_bytes())
        if (observation['formula_units'] != 4 or observation['n_total_images'] != 9
                or observation['climb'] or observation['pressure_GPa'] != 0
                or observation['fmax_target_eV_A'] != .10
                or (not initial and observation['replayed_fmax_eV_A'] > .10)):
            raise ValueError('same ordinary nine-image midpoint target required')
        np.testing.assert_allclose(observation['mechanical_boundary']['reference_cell_A'][:2],
            np.asarray(old_obs['mechanical_boundary']['reference_cell_A'])[:2], atol=1e-12, rtol=0)
        images = []
        for i, row in enumerate(observation['raw_image_evaluations']):
            path = source.parent/f'POSCAR_{i:02d}'
            raw = source.parent/'raw'/f'image_{i:04d}'
            if (sha256(path) != row['snapshot_POSCAR_sha256']
                    or any(row['input_sha256'][n] != h for n, h in CONTRACT.items())
                    or any(sha256(raw/n) != row['input_sha256'][n] for n in ('INPUT', 'KPT', 'STRU'))
                    or sha256(raw/'OUT.ABACUS/running_scf.log') != row['raw_log_sha256']):
                raise ValueError('portable native geometry/physics/output changed')
            parsed = audited_results(raw)
            if any(not np.allclose(parsed[k], v, atol=1e-12, rtol=0) for k, v in (
                    ('energy', row['energy_eV_cell']), ('forces', row['forces_eV_A']),
                    ('stress', row['stress_ASE_voigt_eV_A3']))):
                raise ValueError('portable energy/forces/stress do not replay')
            images.append(read(path, format='vasp'))
        for i in (0, 8):
            old_image = read(repository/RESULTS/channel/f"step_{old_audit['step']:04d}"/f'POSCAR_{i:02d}', format='vasp')
            a, b = observation['raw_image_evaluations'][i], old_obs['raw_image_evaluations'][i]
            if (not same_ordered_geometry(images[i], old_image) or a['input_sha256'] != b['input_sha256']
                    or a['raw_log_sha256'] != b['raw_log_sha256'] or a['energy_eV_cell'] != b['energy_eV_cell']):
                raise ValueError('both starting paths must use the same exact target endpoints')
        transported = original_chart_images(images, orientation, t)
        q, polar, chart, residual = path_coordinates(transported, parent, basis)
        energy = np.array([r['energy_eV_cell'] for r in observation['raw_image_evaluations']])
        relative = (energy-energy[0])*250
        arc = np.asarray(observation['extended_reaction_coordinate_A'])
        registered_peak = audit['sampled_initial_peak_meV_fu'] if initial else audit['sampled_peak_eV_fu']*1000
        if np.any(np.diff(arc) <= 0) or not np.isclose(relative.max(), registered_peak, atol=1e-9):
            raise ValueError('actual energies and generalized arc required')
        records = [dict(channel=channel, initialization='from_1percent', image=i, s=float(arc[i]/arc[-1]),
            relative_energy_meV_fu=float(relative[i]), Q_x_A=float(q[i, 0]), Q_z_A=float(q[i, 2]),
            polar_dz_A=float(polar[i, 2]), open_y_stretch_percent=float(100*(chart['deformation_gradients'][i, 1, 1]-1)),
            residual_atomic_A=float(residual[i])) for i in range(9)]
        old_records = [dict(channel=channel, initialization='from_zero', image=r['image'], s=r['s'],
            relative_energy_meV_fu=r['relative_energy_meV_fu'], Q_x_A=r['Q_x_A'], Q_z_A=r['Q_z_A'],
            polar_dz_A=r['polar_dz_A'], open_y_stretch_percent=100*(r['F_yy']-1), residual_atomic_A=r['residual_atomic_A'])
            for r in old_band['rows']]
        bands.extend([dict(channel=channel, initialization='from_zero', rows=old_records),
                      dict(channel=channel, initialization='from_1percent', rows=records)])
        peak = records[int(relative.argmax())]; old_peak = old_records[old_band['peak_image']]
        outcomes.append(dict(channel=channel, from_zero_peak_meV_fu=old_peak['relative_energy_meV_fu'],
            from_upper_peak_meV_fu=peak['relative_energy_meV_fu'],
            from_zero_Qx_A=old_peak['Q_x_A'], from_upper_Qx_A=peak['Q_x_A'],
            from_zero_open_y_percent=old_peak['open_y_stretch_percent'], from_upper_open_y_percent=peak['open_y_stretch_percent'],
            from_zero_updates=old_audit['step'], from_upper_updates=audit['step'],
            from_zero_fmax_eV_A=old_audit['final_fmax_eV_A'],
            from_upper_fmax_eV_A=audit['replayed_fmax_eV_A'] if initial else audit['final_fmax_eV_A'],
            upper_ordinary_converged=bool(observation['ordinary_residual_pass']),
            upper_terminal_path=not initial,
            upper_SCFs=audit['new_SCFs_in_this_initial_evaluation'] if initial else audit['new_SCF_calls'],
            upper_allocation_core_hours=None if initial else audit['allocation_core_hours'],
            source_observation_sha256=audit['observation_sha256']))
    status = ('initial_seed_comparison_not_terminal_paths' if initial else
              'two_initializations_same_midpoint_compared' if len(channels) == 2 else
              'one_completed_midpoint_channel_compared')
    return dict(status=status,
        bands=bands, outcomes=outcomes, upper_stage='initial_evaluation' if initial else 'ordinary_terminal',
        paired=paired_summary(outcomes) if len(channels) == 2 else None,
        complete_selected_pair_available=len(channels) == 2,
        target_strain=.005, new_DFT_calls=0,
        old_forecast_manifest_sha256=original['forecast_manifest_sha256'],
        not_a_new_independent_prediction_test=True, stationary_TS_certified=False)


def figure(data):
    # Quantitative grid: energies, oxygen rearrangement and cell response form
    # one evidence chain. Each line retains its own actual normalized arc.
    if {b['channel'] for b in data['bands']} != set(JOBS):
        raise ValueError('two-channel figure requires both completed selected channels')
    plt = style()
    from matplotlib.lines import Line2D
    fig, axes = plt.subplots(3, 2, figsize=(183/25.4, 190/25.4))
    fig.subplots_adjust(left=.12, right=.98, bottom=.09, top=.88, wspace=.32, hspace=.36)
    for col, channel in enumerate(('preserving', 'T')):
        for initializer, color, marker in (('from_zero', '#306A94', 'o'), ('from_1percent', '#D48345', 's')):
            rows = next(b['rows'] for b in data['bands'] if b['channel'] == channel and b['initialization'] == initializer)
            s = [r['s'] for r in rows]; kw = dict(color=color, marker=marker, markersize=3.5)
            axes[0, col].plot(s, [r['relative_energy_meV_fu'] for r in rows], **kw)
            axes[1, col].plot(s, [r['Q_x_A'] for r in rows], **kw)
            axes[1, col].plot(s, [r['Q_z_A'] for r in rows], linestyle='--', **kw)
            axes[2, col].plot(s, [r['open_y_stretch_percent'] for r in rows], **kw)
        for ax in axes[:, col]:
            ax.set_xlim(-.035, 1.035); ax.set_xticks([0, .5, 1])
        axes[0, col].set_ylim(-5, 155); axes[0, col].set_yticks([0, 50, 100, 150])
        axes[1, col].set_ylim(-1.17, 1.05); axes[1, col].set_yticks([-1, 0, 1])
        axes[2, col].set_ylim(-1.6, 2.05); axes[2, col].set_yticks([-1, 0, 1, 2])
        axes[2, col].set_xlabel('Normalized path arc\n'+('Preserving flip' if channel == 'preserving' else 'PO to T'))
    axes[0, 0].set_ylabel('Energy (meV/HfO$_2$)')
    axes[1, 0].set_ylabel('Shuffle amplitude (Å)')
    axes[2, 0].set_ylabel('Open-axis stretch (%)')
    upper_label = '1% → 0.5%: initial' if data['upper_stage'] == 'initial_evaluation' else '1% → 0.5%: converged'
    handles = [Line2D([], [], color=c, marker=m, label=label) for c, m, label in
        (('#306A94', 'o', '0% → 0.5%: converged'), ('#D48345', 's', upper_label))]
    # One two-row legend remains above the panel labels; columns are filled
    # vertically by matplotlib, so initialization labels occupy its first row.
    fig.legend(handles=[handles[0], Line2D([], [], color='#555555', label='$Q_x$'),
                        handles[1], Line2D([], [], color='#555555', linestyle='--', label='$Q_z$')],
               loc='upper center', bbox_to_anchor=(.55, .995), ncol=2)
    # Row-major order gives (a,b), (c,d), (e,f) in reading order.
    return finish(fig, axes)


def draw_result(data, output):
    import matplotlib.pyplot as plt
    fig = figure(data)
    for ext in ('pdf', 'png'):
        fig.savefig(output/f'hfo2_midpoint_seed.{ext}', dpi=300, facecolor='white')
    svg = StringIO(); fig.savefig(svg, format='svg', facecolor='white')
    (output/'hfo2_midpoint_seed.svg').write_text('\n'.join(s.rstrip() for s in svg.getvalue().splitlines())+'\n', encoding='utf-8')
    plt.close(fig)


def write_result(data, output, draw=False):
    if draw and {b['channel'] for b in data['bands']} != set(JOBS):
        raise ValueError('do not fabricate a missing channel for the two-channel figure')
    output.mkdir(parents=True, exist_ok=False)
    with (output/'comparison.json').open('x', encoding='utf-8') as handle:
        json.dump(data, handle, indent=2, allow_nan=False); handle.write('\n')
    for name, rows in (('outcomes', data['outcomes']), ('profiles', [r for b in data['bands'] for r in b['rows']])):
        with (output/f'{name}.csv').open('x', encoding='utf-8', newline='') as handle:
            writer = csv.DictWriter(handle, list(rows[0])); writer.writeheader(); writer.writerows(rows)
    if draw:
        draw_result(data, output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repository', type=Path, default=Path(__file__).resolve().parents[4])
    parser.add_argument('--exported', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--draw', action='store_true')
    parser.add_argument('--initial', action='store_true', help='explicitly compare nonterminal step-zero snapshots')
    parser.add_argument('--channel', choices=('preserving', 'T'), help='process one completed channel without a paired-gap claim')
    args = parser.parse_args()
    if args.channel and args.draw:
        parser.error('the two-channel figure needs both channels; omit --draw for a partial result')
    result = compare(args.repository, args.exported, initial=args.initial,
                     channels=(args.channel,) if args.channel else ('preserving', 'T'))
    write_result(result, args.output, args.draw)
    print(json.dumps(dict(outcomes=result['outcomes'], paired=result['paired']), allow_nan=False))
