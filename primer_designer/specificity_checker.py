import json
import os
from dataclasses import dataclass, asdict
from typing import List, Optional
from primer_designer.blast_engine import run_blast, BlastConfig, BlastHit

@dataclass
class SpecificityReport:
    primer_id: str
    total_hits: int
    on_target_hits: List[BlastHit]
    off_target_hits: List[BlastHit]
    pseudogene_hits: List[BlastHit]
    multi_locus_risk: bool
    max_off_target_identity: float
    specificity_score: float
    risk_level: str
    recommendation: str


def _classify_risk(total_hits: int, max_off_target_identity: float, pseudogene: bool, has_on_target: bool) -> str:
    if pseudogene:
        return "CRITICAL"
    if max_off_target_identity >= 90.0:
        return "CRITICAL"
    if max_off_target_identity >= 85.0 or (total_hits > 3 and not has_on_target):
        return "HIGH"
    if max_off_target_identity >= 80.0 or (total_hits > 1 and not has_on_target):
        return "MEDIUM"
    return "LOW"


def _score_specificity(total_hits: int, max_off_target_identity: float, pseudogene: bool, has_on_target: bool) -> float:
    if pseudogene:
        return 20.0
    penalty = min(total_hits * 10.0, 50.0)
    if has_on_target:
        penalty *= 0.5
    score = 100.0 - penalty - max(0.0, max_off_target_identity - 70.0)
    return max(0.0, min(100.0, score))


def run_primer_specificity(primer_seq: str, primer_id: str, cache_dir: str = "cache") -> SpecificityReport:
    config = BlastConfig(program="blastn", database="nt", identity_threshold=70.0, coverage_threshold=0.0, max_hits=20, e_value_cutoff=1e-3, use_cache=True)
    hits = run_blast(primer_seq, config=config, cache_dir=cache_dir)

    on_target: List[BlastHit] = []
    off_target: List[BlastHit] = []
    pseudo: List[BlastHit] = []
    max_identity = 0.0
    max_off_target_identity = 0.0
    for hit in hits:
        if hit.identity_percent > max_identity:
            max_identity = hit.identity_percent
        if hit.organism and "pseudogene" in hit.organism.lower():
            pseudo.append(hit)
        if hit.identity_percent >= 95.0 and hit.query_coverage >= 80.0:
            on_target.append(hit)
        else:
            off_target.append(hit)
            if hit.identity_percent > max_off_target_identity:
                max_off_target_identity = hit.identity_percent

    total_hits = len(hits)
    has_on_target = len(on_target) > 0
    risk = _classify_risk(total_hits, max_off_target_identity, len(pseudo) > 0, has_on_target)
    score = _score_specificity(total_hits, max_off_target_identity, len(pseudo) > 0, has_on_target)
    recommendation = "ACCEPTED" if risk in ("LOW", "MEDIUM") else "REJECTED"

    return SpecificityReport(
        primer_id=primer_id,
        total_hits=total_hits,
        on_target_hits=on_target,
        off_target_hits=off_target,
        pseudogene_hits=pseudo,
        multi_locus_risk=len(off_target) > 1,
        max_off_target_identity=max_identity,
        specificity_score=score,
        risk_level=risk,
        recommendation=recommendation,
    )
