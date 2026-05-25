from dataclasses import dataclass
from typing import List, Optional, Tuple

@dataclass
class StructureRiskReport:
    primer_id: str
    hairpin_detected: bool
    hairpin_sequence: Optional[str]
    hairpin_binding_positions: Optional[Tuple[int, int]]
    hairpin_delta_g: Optional[float]
    self_dimer_detected: bool
    self_dimer_binding_positions: Optional[List[Tuple[int, int]]]
    self_dimer_delta_g: Optional[float]
    cross_dimer_detected: bool
    cross_dimer_binding_positions: Optional[List[Tuple[int, int]]]
    cross_dimer_delta_g: Optional[float]
    three_prime_binding: bool
    three_prime_complement_seq: Optional[str]
    overall_risk: str
    pass_fail: bool
    failure_explanation: Optional[str]


_comp = str.maketrans('ATGCatgc', 'TACGtacg')

def revcomp(s: str) -> str:
    return s.translate(_comp)[::-1]


def _estimate_delta_g(stem_len: int) -> float:
    # simplified heuristic: more negative for longer stems
    return -2.0 * stem_len


def detect_hairpin(seq: str, min_stem: int = 4, max_stem: int = 8, min_loop: int = 3, max_loop: int = 20) -> Optional[Tuple[int, int, int, float]]:
    s = seq.upper()
    L = len(s)
    for i in range(L):
        for stem in range(max_stem, min_stem - 1, -1):
            for loop in range(min_loop, max_loop + 1):
                left_start = i
                left_end = i + stem
                right_start = left_end + loop
                right_end = right_start + stem
                if right_end > L:
                    continue
                left = s[left_start:left_end]
                right = s[right_start:right_end]
                # check complementarity
                if left == revcomp(right):
                    delta_g = _estimate_delta_g(stem)
                    return (left_start + 1, right_end, stem, delta_g)
    return None


def _best_complementarity(a: str, b: str) -> Tuple[int, Tuple[int, int]]:
    # find longest contiguous complementarity between a and b (b reversed)
    ra = a.upper()
    rb = revcomp(b.upper())
    best = 0
    best_pos = (0, 0)
    # sliding window
    for i in range(len(ra)):
        for j in range(len(rb)):
            k = 0
            while i + k < len(ra) and j + k < len(rb) and ra[i + k] == rb[j + k]:
                k += 1
            if k > best:
                best = k
                best_pos = (i + 1, j + 1)
    return best, best_pos


def analyze_primer(primer_seq: str, partner_seq: Optional[str] = None) -> StructureRiskReport:
    hair = detect_hairpin(primer_seq)
    hairpin_detected = hair is not None
    hairpin_seq = None
    hairpin_pos = None
    hairpin_delta = None
    if hair:
        left, right_end, stem, dg = hair
        hairpin_seq = primer_seq[left - 1:right_end]
        hairpin_pos = (left, right_end)
        hairpin_delta = dg

    # self-dimer
    best_self, pos = _best_complementarity(primer_seq, primer_seq)
    self_dimer_detected = best_self >= 4
    self_dimer_positions = [pos] if self_dimer_detected else None
    self_dimer_delta = _estimate_delta_g(best_self) if self_dimer_detected else None

    # cross-dimer
    cross_dimer_detected = False
    cross_positions = None
    cross_delta = None
    three_prime_binding = False
    three_prime_comp = None
    if partner_seq:
        best_cross, pos_cross = _best_complementarity(primer_seq, partner_seq)
        if best_cross >= 4:
            cross_dimer_detected = True
            cross_positions = [pos_cross]
            cross_delta = _estimate_delta_g(best_cross)
            # check 3' binding: if complementarity includes last 5 bases of either primer
            if pos_cross[0] + best_cross - 1 >= len(primer_seq) - 4:
                three_prime_binding = True
                three_prime_comp = partner_seq[::-1][:best_cross]

    # risk
    overall = 'LOW'
    fail = False
    fail_expl = None
    if hairpin_detected and hairpin_delta is not None and hairpin_delta <= -2.0:
        overall = 'MEDIUM'
    if self_dimer_detected and self_dimer_delta is not None and self_dimer_delta <= -3.0:
        overall = 'MEDIUM'
    if cross_dimer_detected and cross_delta is not None and cross_delta <= -5.0:
        overall = 'HIGH'
    if three_prime_binding:
        overall = 'HIGH'
        fail = True
        fail_expl = "3' complementarity detected >= 3 bp"

    return StructureRiskReport(
        primer_id='PRIMER',
        hairpin_detected=hairpin_detected,
        hairpin_sequence=hairpin_seq,
        hairpin_binding_positions=hairpin_pos,
        hairpin_delta_g=hairpin_delta,
        self_dimer_detected=self_dimer_detected,
        self_dimer_binding_positions=self_dimer_positions,
        self_dimer_delta_g=self_dimer_delta,
        cross_dimer_detected=cross_dimer_detected,
        cross_dimer_binding_positions=cross_positions,
        cross_dimer_delta_g=cross_delta,
        three_prime_binding=three_prime_binding,
        three_prime_complement_seq=three_prime_comp,
        overall_risk=overall,
        pass_fail=not fail,
        failure_explanation=fail_expl,
    )
