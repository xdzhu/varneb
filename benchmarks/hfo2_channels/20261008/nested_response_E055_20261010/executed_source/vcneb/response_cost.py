"""Cost of the full visible response dataset, including failures and reuse.

Timing/hardware and source hashes are caller-recorded evidence, not measured
by this API. Transport wall times multiplied by assigned CPU cores are a
resource-cost proxy, not a CPU profiler or prediction-accuracy certificate.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral, Real
import math


@dataclass(frozen=True)
class ResponseEvaluationCost:
    evaluation_id: str
    raw_audit_sha256: str
    outcome: str
    elapsed_seconds: float | None
    allocated_cpu_cores: int | None

    def __post_init__(self):
        if (not isinstance(self.evaluation_id,str) or not self.evaluation_id.strip()
                or not isinstance(self.raw_audit_sha256,str) or len(self.raw_audit_sha256)!=64
                or any(c not in '0123456789abcdef' for c in self.raw_audit_sha256)
                or self.outcome not in ('completed','failed')):
            raise ValueError("actual evaluation ID, raw audit SHA256 and completed/failed outcome required")
        wall=self.elapsed_seconds
        if wall is not None:
            if isinstance(wall,bool) or not isinstance(wall,Real) or not math.isfinite(wall) or wall<0:
                raise ValueError("nonnegative recorded elapsed seconds, or explicitly unknown, required")
            object.__setattr__(self,'elapsed_seconds',float(wall))
        cores=self.allocated_cpu_cores
        if cores is not None:
            if isinstance(cores,bool) or not isinstance(cores,Integral) or cores<1:
                raise ValueError("positive recorded CPU core allocation, or explicitly unknown, required")
            object.__setattr__(self,'allocated_cpu_cores',int(cores))


def recorded_response_dataset_cost(records, *, required_evaluation_ids=()):
    """Count every provided source call, deduplicate shared exact provenance.

    Extra supplied calls and failed preparations remain included, even when
    only a subset of models/columns is eventually used. Missing required IDs
    or inconsistent duplicate records fail closed. Unknown timing/resources
    remain unknown totals, not zero. Cached model accesses are not new calls:
    point them at the same original evaluation ID instead.
    """
    unique={}
    for r in records:
        if not isinstance(r,ResponseEvaluationCost):
            raise TypeError("recorded response-evaluation costs required")
        if r.evaluation_id in unique and unique[r.evaluation_id]!=r:
            raise ValueError("conflicting timing/resource/provenance for one shared evaluation")
        unique[r.evaluation_id]=r
    required=tuple(required_evaluation_ids)
    if any(not isinstance(s,str) or not s.strip() for s in required):
        raise ValueError("explicit required evaluation IDs expected")
    missing=set(required)-set(unique)
    if missing:
        raise ValueError(f"unrecorded required dataset evaluations: {sorted(missing)}")
    rows=list(unique.values())
    unknown_wall=sorted(r.evaluation_id for r in rows if r.elapsed_seconds is None)
    unknown_cpu=sorted(r.evaluation_id for r in rows if r.allocated_cpu_cores is None)
    known_wall=math.fsum(r.elapsed_seconds for r in rows if r.elapsed_seconds is not None)
    known_core_hours=math.fsum(r.elapsed_seconds*r.allocated_cpu_cores/3600 for r in rows
                             if r.elapsed_seconds is not None and r.allocated_cpu_cores is not None)
    if not math.isfinite(known_wall) or not math.isfinite(known_core_hours):
        raise ValueError("nonfinite resource sum; no saturation or zero replacement")
    return {'format_version':1,'unique_source_DFT_evaluations':len(rows),
            'completed_evaluations':sum(r.outcome=='completed' for r in rows),
            'failed_evaluations':sum(r.outcome=='failed' for r in rows),
            'provided_evaluations_not_in_required_subset':len(set(unique)-set(required)),
            'recorded_wall_seconds_sum':None if unknown_wall else known_wall,
            'recorded_transport_cpu_core_hours_sum':None if unknown_wall or unknown_cpu else known_core_hours,
            'known_wall_seconds_partial_sum':known_wall,'known_transport_core_hours_partial_sum':known_core_hours,
            'unknown_elapsed_evaluation_ids':unknown_wall,'unknown_CPU_allocation_evaluation_ids':unknown_cpu,
            'cost_complete':not unknown_wall and not unknown_cpu,'new_DFT_calls_by_this_analysis':0,
            'scope':'All provided visible dataset preparation calls, not only selected columns or successful models',
            'limitations':['source hashes and caller timing/hardware declarations do not perform an external audit',
                           'transport wall times times allocated cores are a cost proxy, not measured CPU utilization',
                           'GPU/node/queue/earlier endpoint costs are not invented; include them separately if applicable',
                           'cost totals do not establish prediction accuracy, acceleration or material validity']}
