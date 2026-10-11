"""Release only the original held +1% reversing starter after audited M.

Reuse the pinned E100 release/seed checks and E103 accepted-child slot policy.
No old file or live source changes: the policy substitution is explicit and
confined to this process. Shared lock/intents remain in the original E100
namespace, so a previous/ambiguous release can never be repeated.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import shlex

RUNS = Path('/public/home/iai806/abacus/agent-runs')
ROOT = RUNS/'20261012-varneb-clamped-reversing-release-E104-r1'
RELEASER = RUNS/'20261011-varneb-clamped-held-flips-E100-r1/release_existing_once_r2.py'
RELEASER_SHA = '0fa713d09056fe344d76241afd3e9490ea0bbcd91ebe2d19db19d77540ae2e9b'
MANAGER = RUNS/'20261012-varneb-clamped-flip-continue-E103-r1/continue_flips_HF.py'
MANAGER_SHA = '6c64c0155800082e28cbbb280511ac761929bdb6ccc5757fa8d68d156b32eac7'
M_AUDIT = RUNS/'20261011-varneb-clamped-child-audit-E093-r1/28727368/terminal/audit_receipt.json'
ZERO_PRESERVING_SHA = '46145da14f24b2b8bd46aa4656f0982650d6b2dc661282f58d809aab8a6d3202'
TARGET = '28709791'
PARENT = '28709790'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_pinned(path, digest, name):
    if sha(path) != digest:
        raise ValueError('reviewed release or accepted-child policy bytes changed')
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def validate_child(body, manager, admitted):
    fields = body.strip().split('|')
    p = manager.paths(PARENT)
    if (len(fields) != 9 or fields[0] != admitted['job_id']
            or fields[1] not in ('RUNNING', 'COMPLETING', 'COMPLETED')
            or fields[2:6] != ['iai806', '32', 'hfacnormal01', 'hfo2-G2-E103-'+PARENT]
            or fields[6:8] != [str(p['stdout']), str(p['stderr'])]):
        raise ValueError('accepted-child native identity/resource/source differs')
    argv = shlex.split(fields[8])
    exports = [v[len('--export='):] for v in argv if v.startswith('--export=')]
    if not argv or argv[0] != 'sbatch' or argv[-1] != str(manager.SCRIPT) or len(exports) != 1:
        raise ValueError('accepted-child native submitted wrapper differs')
    terms = exports[0].split(',')
    pairs = [v.split('=', 1) for v in terms if '=' in v]
    receipt_path = p['terminal']/'audit_receipt.json'
    if (len({k for k, _ in pairs}) != len(pairs)
            or [v for v in terms if '=' not in v] != ['ALL']
            or manager.sha(receipt_path) != admitted['audit_receipt_sha256']):
        raise ValueError('accepted-child ancestry/export identity differs')
    receipt = json.loads(receipt_path.read_text())
    if dict(pairs) != manager.environment(PARENT, receipt, admitted['audit_receipt_sha256']):
        raise ValueError('accepted-child submitted ancestry/physical environment differs')


def reviewed_control():
    legacy = load_pinned(RELEASER, RELEASER_SHA, 'E104_pinned_E100')
    manager = load_pinned(MANAGER, MANAGER_SHA, 'E104_pinned_E103')
    admitted = manager.accepted(PARENT)
    if (legacy.ROOT/'admission.lock' != manager.LOCK
            or manager.sha(manager.SCRIPT) != manager.SCRIPT_SHA):
        raise ValueError('one original shared admission lock and accepted wrapper required')
    native = legacy.run(['sacct', '-X', '-n', '-P', '-j', admitted['job_id'],
        '--format=JobIDRaw,State,User,AllocCPUS,Partition,JobName%100,StdOut%1500,StdErr%1500,SubmitLine%4000'])
    validate_child(native, manager, admitted)
    # PARENT is the original terminal preserving starter, not the live child
    # or the reversing target. Therefore the post-release check can count R.
    legacy.committed_slots = lambda body, require_spare: manager.slots(body, PARENT, spare=require_spare)
    return legacy, admitted, native


def preflight():
    legacy, admitted, native = reviewed_control()
    proof = legacy.preflight(TARGET, require_spare=False)
    return dict(status='read_only_original_reversing_readiness_not_release_or_convergence',
        check_UTC=datetime.now(timezone.utc).isoformat(), target_job=TARGET,
        original_held_seed_proof=proof, accepted_preserving_child=admitted,
        accepted_child_native_source_row=native, executed_adapter_sha256=sha(Path(__file__)),
        original_release_helper_sha256=RELEASER_SHA, accepted_slot_policy_sha256=MANAGER_SHA,
        new_DFT_calls=0, scheduler_writes=0, old_files_or_live_source_changed=False)


def release_once(M_audit_sha256):
    output = ROOT/'release_driver_receipt.json'
    if output.exists():
        raise FileExistsError('existing adapter receipt; inspect without release replay')
    legacy, admitted, native = reviewed_control()
    # A live row, successful cap or timeout cannot replace this physical gate.
    legacy.verified_terminal(M_AUDIT, M_audit_sha256, '28727368')
    result = legacy.release_once(TARGET, ZERO_PRESERVING_SHA)
    proof = dict(status='original_reversing_once_released_with_reviewed_accepted_child_policy',
        target_job=TARGET, check_UTC=datetime.now(timezone.utc).isoformat(),
        M_terminal_audit_sha256=M_audit_sha256, accepted_preserving_child=admitted,
        accepted_child_native_source_row=native, original_release_result=result,
        executed_adapter_sha256=sha(Path(__file__)), original_release_helper_sha256=RELEASER_SHA,
        accepted_slot_policy_sha256=MANAGER_SHA, slot_policy_process_substitution_explicit=True,
        old_files_or_live_source_changed=False, new_submissions=0, dependency_writes=0,
        release_writes=1, new_independent_chains=0, admission_DFT_calls=0)
    with output.open('x') as stream:
        json.dump(proof, stream, indent=2, allow_nan=False)
    return proof


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--release', action='store_true')
    parser.add_argument('--M-audit-sha256')
    args = parser.parse_args()
    if args.release and not args.M_audit_sha256:
        parser.error('release requires the reviewed actual ordinary-converged M audit SHA256')
    print(json.dumps(release_once(args.M_audit_sha256) if args.release else preflight(), indent=2))
