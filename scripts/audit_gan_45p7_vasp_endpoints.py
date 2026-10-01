"""Audit the *production-chain* GaN VASP endpoints against their raw statics.

The signed basin-return runs are separate diagnostics and must not be used as
surrogates for images 0 and 28 of the published 29-image chain.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa

from scripts.audit_gan_ts_basin_followups import residual_stress_kbar


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = ROOT / "paper/VARNEB_CPC/evidence/gan_vasp_endpoint_statics_20260930"
DEFAULT_CHAIN = ROOT / (
    "paper/VARNEB_CPC/evidence/gan_45p7_final_chains_20260927/"
    "gan_vasp_45p7_final_chain.traj"
)
PRESSURE_GPA = 45.7
FORCE_GATE_EV_PER_A = 0.02
STRESS_GATE_KBAR = 2.0
EXPECTED_CHAIN_SHA256 = (
    "952298c1830b293690b8fc4722649147cba236e50e149d07dae45ff75fc2bb7e"
)
EXPECTED_SOURCE_SHA256 = {
    "initial_00_INCAR": "83ba34d4b14ff4ea641bd74dec26e462ea1cc24b575db79a07138ccdfd552772",
    "initial_00_KPOINTS": "b5215f3e7608c27f49d45e8129d3edca2cfaa3ba88b8c389d4b52604715712ee",
    "initial_00_POSCAR": "a59eb4fa87d0dea12cdb5f4856d65f3e4ecce11ffae2add4db2d22e0d10aeb47",
    "initial_00_OUTCAR": "3c6795acf28dad6d00dd33a6f99cd755bf0b0cf8c9c4490bcc1fb3de3005dc87",
    "final_28_POSCAR": "05467c878544f68ca3d34f6a181aa45d8708040451ca09daf1f3462ebcdf60ce",
    "final_28_OUTCAR": "5b69743b754b170807ada26d9ed147de3deb83fc6bc7b8749434b2f7cf98a671",
}
# POTCAR is deliberately not redistributed. These two hashes were checked
# read-only on hf against the production image-15 POTCAR on 2026-09-30.
REMOTE_POTCAR_SHA256 = "f94781ce6cf9b9454c093353320c5e80a2a0aceea6fb8353b33b01ac37b95168"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _input_contract(incar: Path, kpoints: Path) -> dict:
    settings = {}
    for line in incar.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            settings[key.strip().upper()] = value.strip().split("#", 1)[0].strip()
    expected = {
        "ENCUT": 600.0, "EDIFF": 1e-7, "SYMPREC": 1e-4,
        "ISYM": -1, "IBRION": -1, "NSW": 0,
    }
    if any(not np.isclose(float(settings[key]), value, atol=1e-12, rtol=0)
           for key, value in expected.items()):
        raise ValueError("archived endpoint INCAR is not the 600-eV static contract")
    lines = [line.strip() for line in kpoints.read_text(encoding="utf-8").splitlines()]
    if len(lines) < 5 or lines[2].lower() != "gamma" or lines[3] != "8 8 6":
        raise ValueError("archived endpoint KPOINTS is not Gamma 8x8x6")
    return {"ENCUT_eV": 600, "EDIFF_eV": 1e-7, "SYMPREC": 1e-4,
            "ISYM": -1, "k_mesh": [8, 8, 6], "POTCAR_sha256_remote": REMOTE_POTCAR_SHA256}


def audit(source_dir: Path = DEFAULT_SOURCE, chain_path: Path = DEFAULT_CHAIN) -> dict:
    source_dir = Path(source_dir)
    chain_path = Path(chain_path)
    if sha256(chain_path) != EXPECTED_CHAIN_SHA256:
        raise ValueError("VASP final-chain trajectory hash changed")
    for name, expected in EXPECTED_SOURCE_SHA256.items():
        if sha256(source_dir / name) != expected:
            raise ValueError(f"raw endpoint source hash changed: {name}")
    contract = _input_contract(source_dir / "initial_00_INCAR",
                               source_dir / "initial_00_KPOINTS")
    chain = read(chain_path, index=":")
    if len(chain) != 29:
        raise ValueError("expected the archived 29-image VASP chain")

    results = []
    for index, phase, prefix in ((0, "B4", "initial_00"), (28, "B1", "final_28")):
        outcar = source_dir / f"{prefix}_OUTCAR"
        raw_text = outcar.read_text(encoding="utf-8", errors="replace")
        if ("aborting loop because EDIFF is reached" not in raw_text
                or "General timing and accounting informations for this job"
                not in raw_text):
            raise ValueError(f"{phase} raw static lacks SCF or completed-run marker")
        static = read(outcar)
        poscar = read(source_dir / f"{prefix}_POSCAR", format="vasp")
        cached = chain[index]
        if (static.get_chemical_symbols() != cached.get_chemical_symbols()
                or poscar.get_chemical_symbols() != cached.get_chemical_symbols()):
            raise ValueError(f"{phase} atom ordering differs from the chain")
        poscar_cell_diff = float(np.max(np.abs(poscar.cell.array - cached.cell.array)))
        poscar_position_diff = float(np.max(np.abs(poscar.positions - cached.positions)))
        raw_cell_diff = float(np.max(np.abs(static.cell.array - cached.cell.array)))
        raw_position_diff = float(np.max(np.abs(static.positions - cached.positions)))
        energy_diff = float(static.get_potential_energy() - cached.get_potential_energy())
        force_diff = float(np.max(np.abs(static.get_forces() - cached.get_forces())))
        stress_diff = float(np.max(np.abs(static.get_stress(voigt=False)
                                          - cached.get_stress(voigt=False))))
        if (poscar_cell_diff > 1e-10 or poscar_position_diff > 1e-10
                or raw_cell_diff > 1e-8 or raw_position_diff > 1e-5
                or abs(energy_diff) > 1e-8 or force_diff > 1e-5
                or stress_diff > 1e-6):
            raise ValueError(f"{phase} raw static does not match chain endpoint")
        fmax = float(np.linalg.norm(static.get_forces(), axis=1).max())
        stress_kbar = residual_stress_kbar(static.get_stress(voigt=False),
                                           PRESSURE_GPA)
        results.append({
            "phase": phase, "image_index": index,
            "energy_eV_per_cell": float(static.get_potential_energy()),
            "volume_A3": float(static.get_volume()),
            "enthalpy_eV_per_cell": float(static.get_potential_energy()
                                          + PRESSURE_GPA * GPa * static.get_volume()),
            "max_atomic_force_eV_per_A": fmax,
            "max_raw_stress_residual_kbar": stress_kbar,
            "force_pass_0p02eV_per_A": fmax <= FORCE_GATE_EV_PER_A,
            "stress_pass_2kbar": stress_kbar <= STRESS_GATE_KBAR,
            "max_POSCAR_vs_chain_cell_difference_A": poscar_cell_diff,
            "max_POSCAR_vs_chain_position_difference_A": poscar_position_diff,
            "raw_minus_chain_energy_eV": energy_diff,
            "max_raw_minus_chain_force_eV_per_A": force_diff,
            "max_raw_minus_chain_stress_eV_per_A3": stress_diff,
        })
    return {
        "status": "GaN_45p7_VASP_production_endpoint_raw_statics_audited",
        "pressure_GPa": PRESSURE_GPA,
        "total_images": len(chain),
        "contract": contract,
        "endpoints": results,
        "source_sha256": {"trajectory": EXPECTED_CHAIN_SHA256,
                          **EXPECTED_SOURCE_SHA256,
                          "auditor": sha256(Path(__file__))},
        "limitations": [
            "The recorded POTCAR hash was checked on hf but POTCAR is not redistributed or locally reverified.",
            "The 2-kbar stress threshold is an analysis gate, not a VASP input or a VCNEB convergence threshold.",
            "These two original endpoint statics are not the later signed basin-return relaxations.",
            "A raw static residual does not predict the geometry or barrier change after any reoptimization.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--chain", type=Path, default=DEFAULT_CHAIN)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = audit(args.source_dir, args.chain)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "endpoints": report["endpoints"]}))


if __name__ == "__main__":
    main()
