"""Plot the audited GaN B4→B1 endpoint-Γ subspace evolution.

Figure contract: the atomic displacement along the 45.7 GPa path resolves
into endpoint-dependent Γ subspaces while the cell deforms. The plotted Γ
amplitudes are not energy contributions or transition-state eigenmodes.
Archetype: quantitative grid; Python/Matplotlib; 180 × 109 mm; editable SVG,
PDF, PNG/TIFF previews, CSV source data and a machine-readable QA record.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from scripts.analyze_gan_gamma_path_subspaces import audit


plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans']
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['pdf.fonttype'] = 42

COLORS = ('#1C5D99', '#168A85', '#9D5C82')
STRAIN_COLORS = ('#1C5D99', '#BB6A31', '#7F6AB0')


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_source_rows(source_csv: Path, regenerated: list[dict]) -> None:
    """Refuse to draw if the archived CSV differs from raw-mode reconstruction."""

    with source_csv.open(newline='', encoding='utf-8') as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 58 or len(regenerated) != 58:
        raise ValueError('GaN source data must have 29 images for each endpoint basis')
    for existing, fresh in zip(rows, regenerated):
        if set(existing) != set(fresh) or existing['reference_phase'] != fresh['reference_phase']:
            raise ValueError('GaN source data schema or phase differs')
        for key, value in fresh.items():
            if key == 'reference_phase':
                continue
            if not np.isclose(float(existing[key]), float(value), atol=1e-10, rtol=0):
                raise ValueError(f'GaN source data changed at {existing["reference_phase"]} '
                                 f'image {existing["image_index"]}: {key}')


def _column(rows: list[dict], name: str) -> np.ndarray:
    return np.array([float(row[name]) for row in rows], dtype=float)


def plot(directory: Path, path_csv: Path, source_csv: Path,
         archived_audit: Path, output_prefix: Path) -> dict:
    expected = {suffix: output_prefix.with_suffix(suffix)
                for suffix in ('.svg', '.pdf', '.png', '.tiff')}
    qa_path = output_prefix.parent / (output_prefix.name + '_qa.json')
    source_output = output_prefix.parent / (output_prefix.name + '_source_data.csv')
    if any(path.exists() for path in (*expected.values(), qa_path, source_output)):
        raise FileExistsError('refusing to overwrite an existing GaN mode figure or data')
    report, regenerated = audit(directory, path_csv)
    if json.loads(archived_audit.read_text(encoding='utf-8')) != report:
        raise ValueError('archived GaN Gamma audit differs from raw reconstruction')
    validate_source_rows(source_csv, regenerated)
    phases = {phase: [row for row in regenerated if row['reference_phase'] == phase]
              for phase in ('B4', 'B1')}
    b4, b1 = phases['B4'], phases['B1']
    x = _column(b4, 'reaction_coordinate_normalized')
    if (len(x) != 29 or not np.all(np.diff(x) > 0)
            or not np.allclose(x, _column(b1, 'reaction_coordinate_normalized'), atol=1e-12)
            or not np.allclose(_column(b4, 'relative_enthalpy_eV_per_GaN'),
                               _column(b1, 'relative_enthalpy_eV_per_GaN'), atol=1e-12)):
        raise ValueError('endpoint projections do not share one ordered path')

    plt.rcParams.update({
        'font.size': 9.5, 'axes.labelsize': 10, 'xtick.labelsize': 9,
        'ytick.labelsize': 9, 'legend.fontsize': 8.1,
        'axes.linewidth': 0.8, 'axes.spines.top': True,
        'axes.spines.right': True, 'xtick.direction': 'in',
        'ytick.direction': 'in',
    })
    fig, axes = plt.subplots(2, 2, figsize=(7.09, 4.29), sharex=True,
                             gridspec_kw={'hspace': 0.36, 'wspace': 0.29})
    ax_e, ax_b4, ax_b1, ax_strain = axes.flat
    energy = _column(b4, 'relative_enthalpy_eV_per_GaN')
    peak = int(np.argmax(energy))
    ax_e.plot(x, energy, '-', color='#27364A', lw=1.8, zorder=2)
    ax_e.scatter(x, energy, s=13, facecolor='white', edgecolor='#27364A',
                 linewidth=0.8, zorder=3)
    ax_e.scatter([x[peak]], [energy[peak]], marker='D', s=42,
                 facecolor='#C47745', edgecolor='white', linewidth=0.8,
                 zorder=4, label='highest image')
    ax_e.axhline(0, color='#9BA7B0', lw=0.7, ls='--', zorder=0)
    ax_e.set_ylabel(r'$\Delta H$ (eV/GaN)')
    ax_e.legend(loc='upper left', frameon=True, framealpha=0.85,
                facecolor='white', edgecolor='#B9C4CD')

    for phase, ax, rows in (('B4', ax_b4, b4), ('B1', ax_b1, b1)):
        groups = report['phases'][phase]['groups_ranked_by_path_amplitude']
        for rank, (group, color) in enumerate(zip(groups, COLORS), start=1):
            q = _column(rows, f'group{rank}_Q_norm_d0p01_sqrt_amu_A')
            modes = ','.join(str(i + 1) for i in group['mode_indices'])
            label = rf'$\Gamma_{{{modes}}}$, {group["frequency_cm1_at_0p01"]:.0f} cm$^{{-1}}$'
            ax.plot(x, q, color=color, lw=1.6, marker='o', ms=2.5,
                    markevery=4, label=label)
        ax.set_ylabel(fr'$|Q|_{{{phase}\,\Gamma}}$ ($\sqrt{{\mathrm{{amu}}}}\,\AA$)')
        ax.legend(loc='upper left' if phase == 'B4' else 'upper right',
                  frameon=True, framealpha=0.85, facecolor='white',
                  edgecolor='#B9C4CD')

    for key, color, label, style in (
        ('symmetric_strain_xx', STRAIN_COLORS[0], r'$\eta_{xx}$', '-'),
        ('symmetric_strain_yy', STRAIN_COLORS[1], r'$\eta_{yy}$', '--'),
        ('symmetric_strain_zz', STRAIN_COLORS[2], r'$\eta_{zz}$', '-.'),
    ):
        ax_strain.plot(x, 100 * _column(b4, key), color=color,
                       lw=1.7, ls=style, label=label)
    ax_strain.axhline(0, color='#9BA7B0', lw=0.7, zorder=0)
    ax_strain.set_ylabel('Strain vs B4 (%)')
    ax_strain.legend(loc='upper left', ncol=3, frameon=True, framealpha=0.85,
                     facecolor='white', edgecolor='#B9C4CD', columnspacing=0.7)

    for index, ax in enumerate(axes.flat):
        ax.set_xlim(0, 1)
        ax.set_xticks(np.linspace(0, 1, 6))
        ax.tick_params(top=True, right=True, length=3.5, width=0.8)
        ax.text(-0.125, 1.05, f'({chr(97 + index)})', transform=ax.transAxes,
                fontsize=11.5, fontweight='normal', ha='left', va='bottom')
    for ax in axes[1]:
        ax.set_xlabel('Normalized path coordinate')
    fig.subplots_adjust(left=0.095, right=0.975, bottom=0.14, top=0.94)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(expected['.svg'], bbox_inches='tight')
    fig.savefig(expected['.pdf'], bbox_inches='tight')
    fig.savefig(expected['.png'], bbox_inches='tight', dpi=300)
    fig.savefig(expected['.tiff'], bbox_inches='tight', dpi=600)
    plt.close(fig)
    shutil.copyfile(source_csv, source_output)
    qa = {
        'status': 'research_figure_endpoint_Gamma_subspaces_not_TS_certificate',
        'core_conclusion': ('Atomic path changes occupy audited endpoint-Gamma '
                            'subspaces and coexist with substantial cell strain.'),
        'figure_archetype': 'quantitative_grid',
        'backend': 'Python/Matplotlib',
        'intended_size_mm': [180, 109],
        'n_distinct_DFT_images': 29,
        'n_Gamma_real_space_cells': [1, 1, 1],
        'highest_image_index_not_certified_TS': peak,
        'highest_image_enthalpy_eV_per_GaN': float(energy[peak]),
        'panel_map': {'a': 'original path enthalpy at 45.7 GPa',
                      'b': 'B4 endpoint-Gamma optical-subspace amplitudes',
                      'c': 'B1 endpoint-Gamma optical-subspace amplitudes',
                      'd': 'symmetric strain relative to B4'},
        'source_sha256': {'audited_Gamma_report': sha256(archived_audit),
                          'mode_source_csv': sha256(source_csv),
                          'path_source_csv': sha256(path_csv),
                          'script': sha256(Path(__file__))},
        'exports_sha256': {path.name: sha256(path) for path in expected.values()},
        'source_data_sha256': sha256(source_output),
        'review_risks': [
            'Endpoint Gamma modes are not full-variable-cell TS modes.',
            'Subspace amplitudes are not mode-resolved energy contributions.',
            'Frequencies and amplitudes use the original 600 eV path; do not mix with 1000 eV TS diagnostic enthalpies.',
            'Path highest image is a candidate, not a certified first-order saddle.',
        ],
    }
    qa_path.write_text(json.dumps(qa, indent=2) + '\n', encoding='utf-8')
    return qa


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('directory', 'path-csv', 'source-csv', 'archived-audit', 'output-prefix'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(plot(args.directory, args.path_csv, args.source_csv,
                          args.archived_audit, args.output_prefix)))


if __name__ == '__main__':
    main()
