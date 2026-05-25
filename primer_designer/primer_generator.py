from dataclasses import dataclass
from typing import List, Optional, Tuple
from primer_designer.thermodynamics import thermo_profile


@dataclass
class PrimerPair:
    forward: str
    reverse: str
    f_start: int
    f_end: int
    r_start: int
    r_end: int


@dataclass
class PrimerCandidate:
    candidate_id: str
    direction: str
    sequence: str
    length: int
    gc_content: float
    tm_basic: float
    tm_advanced: float
    gc_clamp: bool
    repeat_detected: bool
    repeat_sequence: Optional[str]
    is_accepted: bool
    rejection_reasons: List[str]


def _revcomp(seq: str) -> str:
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "a": "t", "t": "a", "g": "c", "c": "g"}
    return "".join(comp.get(b, b) for b in reversed(seq))


def generate_primer_pair(seq: str, length: int = 20, offset_from_ends: int = 0) -> PrimerPair:
    if len(seq) < length * 2 + offset_from_ends:
        raise ValueError("Sequence too short to generate primer pair with given length and offset")
    f_start = 1 + offset_from_ends
    f_end = f_start + length - 1
    r_end = len(seq) - offset_from_ends
    r_start = r_end - length + 1
    f_seq = seq[f_start - 1:f_end]
    r_seq = seq[r_start - 1:r_end]
    r_revcomp = _revcomp(r_seq)
    return PrimerPair(f_seq, r_revcomp, f_start, f_end, r_start, r_end)


def _gc_content(seq: str) -> float:
    s = seq.upper()
    g = s.count('G')
    c = s.count('C')
    total = sum(s.count(b) for b in 'ATGC')
    return 0.0 if total == 0 else 100.0 * (g + c) / total


def _has_homopolymer(seq: str, max_run: int = 3) -> Tuple[bool, str]:
    # detect runs of same base longer than max_run
    last = ''
    run = 0
    for ch in seq:
        if ch == last:
            run += 1
        else:
            run = 1
            last = ch
        if run > max_run:
            return True, last * run
    return False, None


def _has_dinuc_repeat(seq: str, repeat_len: int = 3) -> Tuple[bool, str]:
    # detect dinucleotide repeats like ATATAT of given repeat_len (number of repeats)
    sequ = seq.upper()
    for i in range(len(sequ) - 2 * repeat_len + 1):
        motif = sequ[i:i+2]
        if motif * repeat_len == sequ[i:i+2*repeat_len]:
            return True, motif * repeat_len
    return False, None


def _gc_clamp(seq: str) -> bool:
    return seq[-1].upper() in ('G', 'C')


def validate_candidate(seq: str, length_min: int = 18, length_max: int = 25, gc_min: float = 40.0, gc_max: float = 60.0) -> Tuple[bool, List[str]]:
    reasons: List[str] = []
    L = len(seq)
    if L < length_min or L > length_max:
        reasons.append(f"Length {L} bp outside allowed range ({length_min}-{length_max})")

    gc = _gc_content(seq)
    if gc < gc_min or gc > gc_max:
        reasons.append(f"GC% {gc:.1f} outside allowed range ({gc_min}-{gc_max})")

    homopoly, homseq = _has_homopolymer(seq, max_run=3)
    if homopoly:
        reasons.append(f"Homopolymer run detected: {homseq}")

    dinuc, dinseq = _has_dinuc_repeat(seq, repeat_len=3)
    if dinuc:
        reasons.append(f"Dinucleotide repeat detected: {dinseq}")

    if not _gc_clamp(seq):
        reasons.append("No GC clamp at 3' end (G/C recommended)")

    # thermodynamics
    tp = thermo_profile(seq)
    # basic Tm used for quick filtering
    if tp.tm_advanced <= 0:
        reasons.append("Invalid Tm calculation")

    is_accepted = len(reasons) == 0
    return is_accepted, reasons


def generate_reverse_candidates_from_region(region_seq: str, length_min: int = 18, length_max: int = 25) -> List[PrimerCandidate]:
    """Generate reverse primer candidates as reverse-complements of 3'-end subsequences.

    The reverse primer binds the antisense strand, so we take subsequences from the
    3' portion of the conserved region and reverse-complement them.
    """
    candidates: List[PrimerCandidate] = []
    cid = 1
    for L in range(length_min, length_max + 1):
        for start in range(0, len(region_seq) - L + 1):
            subseq = region_seq[start:start + L]
            rc_seq = _revcomp(subseq)
            tp = thermo_profile(rc_seq)
            is_ok, reasons = validate_candidate(rc_seq, length_min, length_max)
            pc = PrimerCandidate(
                candidate_id=f"R-{cid:04d}",
                direction="REVERSE",
                sequence=rc_seq,
                length=L,
                gc_content=tp.gc_percent,
                tm_basic=tp.tm_basic,
                tm_advanced=tp.tm_advanced,
                gc_clamp=_gc_clamp(rc_seq),
                repeat_detected=_has_homopolymer(rc_seq)[0] or _has_dinuc_repeat(rc_seq)[0],
                repeat_sequence=_has_homopolymer(rc_seq)[1] or _has_dinuc_repeat(rc_seq)[1],
                is_accepted=is_ok,
                rejection_reasons=reasons,
            )
            candidates.append(pc)
            cid += 1
    return candidates


def generate_candidates_from_region(region_seq: str, length_min: int = 18, length_max: int = 25) -> List[PrimerCandidate]:
    candidates: List[PrimerCandidate] = []
    cid = 1
    for L in range(length_min, length_max + 1):
        for start in range(0, len(region_seq) - L + 1):
            subseq = region_seq[start:start+L]
            tp = thermo_profile(subseq)
            is_ok, reasons = validate_candidate(subseq, length_min, length_max)
            pc = PrimerCandidate(
                candidate_id=f'P-{cid:04d}',
                direction='FORWARD',
                sequence=subseq,
                length=L,
                gc_content=tp.gc_percent,
                tm_basic=tp.tm_basic,
                tm_advanced=tp.tm_advanced,
                gc_clamp=_gc_clamp(subseq),
                repeat_detected=_has_homopolymer(subseq)[0] or _has_dinuc_repeat(subseq)[0],
                repeat_sequence=_has_homopolymer(subseq)[1] or _has_dinuc_repeat(subseq)[1],
                is_accepted=is_ok,
                rejection_reasons=reasons,
            )
            candidates.append(pc)
            cid += 1
    return candidates

