import math
from dataclasses import dataclass
from typing import Dict, List, Tuple

from primer_designer.msa_engine import AlignmentMatrix


@dataclass
class ConservedRegion:
    region_id: str
    start_position: int
    end_position: int
    length: int
    conservation_score: float
    gap_frequency: float
    entropy_mean: float       # mean Shannon entropy across columns (bits; 0 = fully conserved)
    consensus_sequence: str
    per_organism_sequences: Dict[str, str]
    is_suitable_forward: bool
    is_suitable_reverse: bool
    evidence_table: str


def _shannon_entropy(counts: Dict[str, int], total_non_gap: int) -> float:
    """Shannon entropy in bits for a single MSA column.

    Returns 0.0 for perfectly conserved columns (one base dominates);
    returns log2(4) ≈ 2.0 for completely random columns.
    """
    if total_non_gap == 0:
        return 0.0
    h = 0.0
    for count in counts.values():
        if count > 0:
            p = count / total_non_gap
            h -= p * math.log2(p)
    return h


def _column_stats(aligned_seqs: List[str], col: int) -> Tuple[float, float, str, Dict[str, int], float]:
    """Return (conservation, gap_freq, consensus_base, base_counts, entropy) for one column."""
    counts: Dict[str, int] = {}
    gap = 0
    total = len(aligned_seqs)

    for s in aligned_seqs:
        if col >= len(s):
            gap += 1
            continue
        ch = s[col].upper()
        if ch in ("-", "."):
            gap += 1
        else:
            counts[ch] = counts.get(ch, 0) + 1

    non_gap = total - gap
    most_common = 0
    consensus = "N"
    for k, v in counts.items():
        if v > most_common:
            most_common = v
            consensus = k

    conservation = (most_common / non_gap) if non_gap > 0 else 0.0
    gap_freq = gap / total
    entropy = _shannon_entropy(counts, non_gap)
    return conservation, gap_freq, consensus, counts, entropy


def find_conserved_regions(
    aln: AlignmentMatrix,
    min_conservation: float = 0.9,
    max_gap: float = 0.05,
    min_length: int = 18,
) -> List[ConservedRegion]:
    """Detect conserved regions in a multiple sequence alignment.

    A column is included when conservation >= min_conservation AND
    gap_frequency <= max_gap. Contiguous runs of qualifying columns that
    span at least min_length positions become ConservedRegion objects.

    The returned list is ordered by start_position (left to right in the alignment).
    """
    seqs = aln.aligned_sequences
    if not seqs:
        return []

    L = aln.alignment_length or (len(seqs[0]) if seqs else 0)

    # Per-column statistics
    cols = []
    for i in range(L):
        cons, gapf, cons_base, counts, entropy = _column_stats(seqs, i)
        cols.append({
            "pos": i + 1,
            "conservation": cons,
            "gap": gapf,
            "consensus": cons_base,
            "counts": counts,
            "entropy": entropy,
        })

    regions: List[ConservedRegion] = []
    current = None
    rid = 1

    def _emit_region(start: int, end: int) -> ConservedRegion:
        nonlocal rid
        span = range(start - 1, end)  # 0-based slice
        cons_scores = [cols[j]["conservation"] for j in span]
        gap_scores = [cols[j]["gap"] for j in span]
        entropy_vals = [cols[j]["entropy"] for j in span]
        consensus_seq = "".join(cols[j]["consensus"] for j in span)
        per_org = {f"seq_{idx + 1}": seq[start - 1:end] for idx, seq in enumerate(seqs)}
        evidence_lines = [f"seq_{idx + 1}: {seq[start - 1:end]}" for idx, seq in enumerate(seqs)]

        region = ConservedRegion(
            region_id=f"CR-{rid:03d}",
            start_position=start,
            end_position=end,
            length=end - start + 1,
            conservation_score=100.0 * (sum(cons_scores) / len(cons_scores)),
            gap_frequency=100.0 * (sum(gap_scores) / len(gap_scores)),
            entropy_mean=sum(entropy_vals) / len(entropy_vals),
            consensus_sequence=consensus_seq,
            per_organism_sequences=per_org,
            is_suitable_forward=True,
            is_suitable_reverse=True,
            evidence_table="\n".join(evidence_lines),
        )
        rid += 1
        return region

    for i, c in enumerate(cols):
        ok = (c["conservation"] >= min_conservation) and (c["gap"] <= max_gap)
        if ok:
            if current is None:
                current = {"start": i + 1, "end": i + 1}
            else:
                current["end"] = i + 1
        else:
            if current is not None:
                length = current["end"] - current["start"] + 1
                if length >= min_length:
                    regions.append(_emit_region(current["start"], current["end"]))
                current = None

    # Handle trailing conserved run
    if current is not None:
        length = current["end"] - current["start"] + 1
        if length >= min_length:
            regions.append(_emit_region(current["start"], current["end"]))

    return regions
