"""Material evidence update, zero DFT. Reuse existing audited mode analysis."""
import argparse
import json
from pathlib import Path

from scripts.analyze_hfo2_network_update import analyze_case
from scripts.analyze_hfo2_channel_network import read_evaluated_observation
from scripts.audit_hfo2_static_replica import sha256, audited_results
from vcneb.polarization import sampled_band_gap


def audit(specification, *, raw_on_HF=False):
    result = analyze_case(specification)
    if len(result["channels"]) != 4 or not all(c["ordinary_residual_passed"] for c in result["channels"]):
        raise ValueError("four actual ordinary terminal passes required")
    # analyze_case's historical last limitation described its original stage.
    # Replace only that stage-specific statement; do not alter its raw methods.
    result["limitations"][-1] = "four ordinary residual passes; full G1, branch stability, TS and strain prediction gates pending"
    result["terminal_update_driver_sha256"] = sha256(Path(__file__))
    result["raw_on_HF"] = raw_on_HF
    if raw_on_HF:
        spec = json.loads(specification.read_text())
        for entry, channel in zip(spec["channels"], result["channels"]):
            images, record, digest = read_evaluated_observation(specification.parent / entry["observation"])
            bands = []
            for i, (atoms, source) in enumerate(zip(images, record["raw_image_evaluations"])):
                directory = Path(source["raw_source"])
                raw = audited_results(directory)
                if abs(raw["energy"] - atoms.get_potential_energy()) > 1e-10:
                    raise ValueError("raw terminal energy changed")
                log = directory / "OUT.ABACUS/running_scf.log"
                if sha256(log) != source["raw_log_sha256"]:
                    raise ValueError("raw terminal log changed")
                p = directory / "OUT.ABACUS/istate.info"
                bands.append({"source_image_index": i, "raw_source": str(directory),
                              "istate_sha256": sha256(p), **sampled_band_gap(p.read_text(), 48)})
            channel["sampled_gap_audit"] = bands
            channel["minimum_sampled_gap_eV"] = min(b["sampled_indirect_gap_eV"] for b in bands)
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--specification", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--raw-on-HF", action="store_true")
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError("refusing existing terminal audit")
    result = audit(a.specification, raw_on_HF=a.raw_on_HF)
    a.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps({"records": result["existing_image_records"], "ordinary_passes": 4, "new_DFT": 0,
                      "sampled_gap_minima": {c["name"]: c.get("minimum_sampled_gap_eV") for c in result["channels"]}}))
