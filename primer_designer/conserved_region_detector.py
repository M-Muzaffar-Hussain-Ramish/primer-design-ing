from dataclasses import dataclass
from typing import List, Dict
from primer_designer.msa_engine import AlignmentMatrix

@dataclass
class ConservedRegion:
    region_id: str
    start_position: int
    end_position: int
    length: int
    conservation_score: float
    gap_frequency: float
    entropy_mean: float
    consensus_sequence: str
    per_organism_sequences: Dict[str, str]
    is_suitable_forward: bool
    is_suitable_reverse: bool
    evidence_table: str


def _column_stats(aligned_seqs: List[str], col: int):
    counts = {}
    gap = 0
    total = len(aligned_seqs)
    for s in aligned_seqs:
        ch = s[col]
        if ch == '-' or ch == '.':
            gap += 1
            continue
        counts[ch] = counts.get(ch, 0) + 1
    non_gap = total - gap
    most_common = 0
    consensus = 'N'
    for k, v in counts.items():
        if v > most_common:
            most_common = v
            consensus = k
    conservation = (most_common / non_gap) if non_gap > 0 else 0.0
    gap_freq = gap / total
    return conservation, gap_freq, consensus, counts


def find_conserved_regions(aln: AlignmentMatrix, min_conservation: float = 0.9, max_gap: float = 0.05, min_length: int = 18) -> List[ConservedRegion]:
    seqs = aln.aligned_sequences
    if not seqs:
        return []
    L = len(seqs[0])
    # compute per-column stats
    cols = []
    for i in range(L):
        cons, gapf, cons_base, counts = _column_stats(seqs, i)
        cols.append({'pos': i+1, 'conservation': cons, 'gap': gapf, 'consensus': cons_base, 'counts': counts})

    regions = []
    current = None
    rid = 1
    for i, c in enumerate(cols):
        ok = (c['conservation'] >= min_conservation) and (c['gap'] <= max_gap)
        if ok:
            if current is None:
                current = {'start': i+1, 'end': i+1}
            else:
                current['end'] = i+1
        else:
            if current is not None:
                length = current['end'] - current['start'] + 1
                if length >= min_length:
                    # build region
                    start = current['start']
                    end = current['end']
                    cons_scores = [cols[j-1]['conservation'] for j in range(start, end+1)]
                    gap_scores = [cols[j-1]['gap'] for j in range(start, end+1)]
                    consensus_seq = ''.join(cols[j-1]['consensus'] for j in range(start, end+1))
                    per_org = {f'seq_{idx+1}': seq[start-1:end] for idx, seq in enumerate(seqs)}
                    evidence_lines = []
                    for idx, seq in enumerate(seqs):
                        evidence_lines.append(f'seq_{idx+1}: {seq[start-1:end]}')
                    evidence = '\n'.join(evidence_lines)
                    region = ConservedRegion(
                        region_id=f'CR-{rid:03d}',
                        start_position=start,
                        end_position=end,
                        length=length,
                        conservation_score=100.0 * (sum(cons_scores) / len(cons_scores)),
                        gap_frequency=100.0 * (sum(gap_scores) / len(gap_scores)),
                        entropy_mean=0.0,
                        consensus_sequence=consensus_seq,
                        per_organism_sequences=per_org,
                        is_suitable_forward=True,
                        is_suitable_reverse=True,
                        evidence_table=evidence,
                    )
                    regions.append(region)
                    rid += 1
                current = None
    # check tail
    if current is not None:
        length = current['end'] - current['start'] + 1
        if length >= min_length:
            start = current['start']
            end = current['end']
            cons_scores = [cols[j-1]['conservation'] for j in range(start, end+1)]
            gap_scores = [cols[j-1]['gap'] for j in range(start, end+1)]
            consensus_seq = ''.join(cols[j-1]['consensus'] for j in range(start, end+1))
            per_org = {f'seq_{idx+1}': seq[start-1:end] for idx, seq in enumerate(seqs)}
            evidence_lines = []
            for idx, seq in enumerate(seqs):
                evidence_lines.append(f'seq_{idx+1}: {seq[start-1:end]}')
            evidence = '\n'.join(evidence_lines)
            region = ConservedRegion(
                region_id=f'CR-{rid:03d}',
                start_position=start,
                end_position=end,
                length=length,
                conservation_score=100.0 * (sum(cons_scores) / len(cons_scores)),
                gap_frequency=100.0 * (sum(gap_scores) / len(gap_scores)),
                entropy_mean=0.0,
                consensus_sequence=consensus_seq,
                per_organism_sequences=per_org,
                is_suitable_forward=True,
                is_suitable_reverse=True,
                evidence_table=evidence,
            )
            regions.append(region)
    return regions
