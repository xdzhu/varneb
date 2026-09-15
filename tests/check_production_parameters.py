"""Guard the material-specific ABACUS production parameter contract."""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    bto = (ROOT / "cluster" / "hf_batio3_vcneb_parallel.slurm").read_text(encoding="utf-8")
    hfo2 = (ROOT / "cluster" / "hf_hfo2_vcneb_distributed.slurm").read_text(encoding="utf-8")
    for label, text, species in (("BTO", bto, ("Ba", "Ti", "O")), ("HfO2", hfo2, ("Hf", "O"))):
        if "ecutwfc=${ECUTWFC:-100}" not in text:
            raise SystemExit(f"{label} template does not default to 100 Ry")
        if "Orb-DZP-10au" not in text:
            raise SystemExit(f"{label} template does not select Orb-DZP-10au")
        for element in species:
            if f"{element}_gga_10au_100Ry" not in text:
                raise SystemExit(f"{label} template lacks the 10-au {element} orbital")

    for fragment in (
        "mode_cell_scale=${MODE_CELL_SCALE:-}",
        'mode_args+=(--mode-cell-scale "${mode_cell_scale}")',
    ):
        if fragment not in hfo2:
            raise SystemExit("HfO2 distributed template does not forward MODE_CELL_SCALE")

    smoke = (ROOT / "examples" / "run_abacus_single_image_smoke.py").read_text(encoding="utf-8")
    endpoint = (ROOT / "examples" / "relax_abacus_endpoint.py").read_text(encoding="utf-8")
    for label, text in (("HfO2 smoke", smoke), ("HfO2 endpoint fallback", endpoint)):
        if "default=100.0" not in text or "Orb-DZP-10au" not in text:
            raise SystemExit(f"{label} defaults drifted from the 100-Ry/10-au policy")
        if "default=1e-8" not in text or (
            "default=[2, 2, 2]" not in text and '"kpts": [2, 2, 2]' not in text
        ):
            raise SystemExit(f"{label} defaults drifted from the production SCF/k-point policy")

    production_readme = (ROOT / "README_VCNEB.md").read_text(encoding="utf-8")
    if "--ecutwfc 100 --kpts 2 2 2 --scf-thr 1e-8" not in production_readme:
        raise SystemExit("README production HfO2 command is missing the high-accuracy settings")
    print("production_parameters_regression=ok")


if __name__ == "__main__":
    main()
