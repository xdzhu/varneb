"""Actual seen M strains, using the same guarded pair builder and renderer."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from scripts import plot_hfo2_clamped_comparison as paired
from scripts.plot_hfo2_clamped_M import (
    BASE, MODE_REPORT, MODE_CANONICAL_SHA, OBSERVATION as ZERO_OBSERVATION,
    OBSERVATION_SHA as ZERO_OBSERVATION_SHA, AUDIT as ZERO_AUDIT,
    AUDIT_SHA as ZERO_AUDIT_SHA, sha256,
)
from scripts.plot_hfo2_clamped_T_strain import build_strain_data

REPORT = BASE/'clamped_M_strain_reference_E105_20261012/displacement.json'
GRADIENT = REPORT.parent/'gradient.json'
TERMINAL = BASE/'clamped_reversing_release_E104_20261012/plus1_M_terminal'
OBSERVATION = TERMINAL/'observations/step_0022/observation.json'
AUDIT = TERMINAL/'audit_receipt.json'
REPORT_SHA = '244055c4b1fcb9c0c810db0ebfb38a78fe385e60f6ef0dfa5a8cdde4b73bca5e'
GRADIENT_SHA = 'fa79242ad672d165a96dfce09d0603e27d5ea2817812472f62c7a2379d6ca3d1'
OBSERVATION_SHA = '634de4416788c90bfe9632af7a9fd1c7dc9a0da5bdb6c9c2f1960a277994277c'
AUDIT_SHA = 'b86cdb72eab9e3d8ff2125e16e6fbb2d0954b34d54c0f05099a581c066108389'


def build_data(repository):
    specs = [
        (0., MODE_REPORT, paired.M_GRADIENT, ZERO_OBSERVATION, ZERO_AUDIT,
         MODE_CANONICAL_SHA, paired.M_GRADIENT_SHA, ZERO_OBSERVATION_SHA,
         ZERO_AUDIT_SHA, '28722320', 13),
        (.01, REPORT, GRADIENT, OBSERVATION, AUDIT, REPORT_SHA, GRADIENT_SHA,
         OBSERVATION_SHA, AUDIT_SHA, '28727368', 22),
    ]
    data = build_strain_data(repository, specs, channel='M')
    data['energy_limits_meV_fu'] = [-135., 125.]
    if any(not -135. < row['energy_relative_PO_meV_fu'] < 125.
           for band in data['bands'] for row in band['rows']):
        raise ValueError('declared paired energy limits do not enclose the native samples')
    data['strain_pair_builder_sha256'] = data['comparison_script_sha256']
    data['comparison_script_sha256'] = sha256(Path(__file__))
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repository', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = paired.export(build_data(args.repository), args.output, stem_name='hfo2_clamped_M_strain')
    print(json.dumps(dict(records=result['image_records'], grid_aligned=result['grid_aligned'],
                          new_DFT_calls=0, independent_prediction=False)))


if __name__ == '__main__':
    main()
