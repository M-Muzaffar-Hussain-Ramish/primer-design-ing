from dataclasses import dataclass
from typing import Optional

@dataclass
class ScoredPrimer:
    primer_id: str
    conservation_score: float
    thermodynamic_score: float
    specificity_score: float
    secondary_structure_score: float
    sequence_quality_score: float
    composite_score: float
    confidence_rating: str
    risk_level: str


def _normalize(value: float, low: float = 0.0, high: float = 100.0) -> float:
    if value <= low:
        return 0.0
    if value >= high:
        return 100.0
    return 100.0 * (value - low) / (high - low)


def score_primer(
    primer_id: str,
    conservation_score: float,
    tm_basic: float,
    tm_advanced: float,
    specificity_score: float,
    secondary_structure_risk: str,
    gc_content: float,
    length: int,
) -> ScoredPrimer:
    # Conservation: already 0-100
    conservation = max(0.0, min(100.0, conservation_score))

    # Thermodynamic score: ideal Tm 55-65 and Tm diff within 2
    tm_mean = (tm_basic + tm_advanced) / 2.0
    tm_bonus = max(0.0, 100.0 - abs(tm_mean - 60.0) * 3.0)
    tm_score = _normalize(tm_bonus, 0.0, 100.0)

    # Specificity score used as raw 0-100
    specificity = max(0.0, min(100.0, specificity_score))

    # Secondary structure: convert risk label to score
    ss_score = {
        'LOW': 100.0,
        'MEDIUM': 60.0,
        'HIGH': 30.0,
        'FAIL': 0.0,
    }.get(secondary_structure_risk.upper(), 50.0)

    # Sequence quality: based on GC content and length
    gc_penalty = max(0.0, abs(gc_content - 50.0) * 2.0)
    length_penalty = 0.0 if 20 <= length <= 22 else 10.0
    quality = max(0.0, 100.0 - gc_penalty - length_penalty)

    # Composite weights
    composite = (
        conservation * 0.25
        + tm_score * 0.20
        + specificity * 0.25
        + ss_score * 0.15
        + quality * 0.15
    )

    if composite >= 90.0:
        confidence = 'HIGH'
        risk = 'LOW'
    elif composite >= 75.0:
        confidence = 'MEDIUM'
        risk = 'MEDIUM'
    else:
        confidence = 'LOW'
        risk = 'HIGH'

    return ScoredPrimer(
        primer_id=primer_id,
        conservation_score=conservation,
        thermodynamic_score=tm_score,
        specificity_score=specificity,
        secondary_structure_score=ss_score,
        sequence_quality_score=quality,
        composite_score=round(composite, 2),
        confidence_rating=confidence,
        risk_level=risk,
    )
