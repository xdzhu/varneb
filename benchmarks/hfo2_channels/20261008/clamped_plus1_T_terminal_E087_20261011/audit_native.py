"""Run the existing pinned terminal auditor for the actual +1% T job only.

No DFT, admission, dependency edit or production-source write. Output must
be fresh. Imported audit modules are pinned separately from the old runtime.
"""
import hashlib
import json
from pathlib import Path
import sys

BASE = Path('/public/home/iai806/abacus/agent-runs')
RUNTIME = BASE/'20261010-varneb-nested-controls-E055-r1/source'
sys.path.insert(0, str(RUNTIME))
import scripts

AUDIT_SOURCE = BASE/'20261011-varneb-clamped-T-terminal-E075-source-r1/scripts'
scripts.__path__ = [str(AUDIT_SOURCE), *list(scripts.__path__)]
import scripts.audit_hfo2_clamped_terminal as auditor
import scripts.export_hfo2_clamped_observation as exporter
import scripts.analyze_hfo2_clamped_residual as residual

PINNED = (
    (auditor, '25079147be017329a435749fe94354c6c5f1740889bf2866c21efd40e8b320c0'),
    (exporter, '83b41f1f1ac8982db40606e279c5fccba3fac44c786c2238a42d5841cd723cea'),
    (residual, 'a2ed0dc5389524b10b77f3329945ab04ad188fe99e57647b574df860a9920d84'),
)
if any(hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest() != digest
       for module, digest in PINNED):
    raise ValueError('actual imported terminal audit modules changed')

if __name__ == '__main__':
    receipt = auditor.audit(
        BASE/'20261010-varneb-clamped-G2-E060-r1/PO_to_T/band', '28709788',
        RUNTIME/'cluster/hf_hfo2_clamped_chain_pilot_20261010.slurm', RUNTIME,
        BASE/'varneb-E055-clean-86d2637.tar',
        'cbe6552a3fa4237fec2098c08c8056735a9a95df4e18691d09fed6b4f24a613c',
        BASE/'20261011-varneb-clamped-plus1-T-terminal-E086-r1/terminal')
    print(json.dumps({k: receipt[k] for k in
                     ('status', 'job_id', 'terminal_step', 'fresh_interior_SCFs')}))
