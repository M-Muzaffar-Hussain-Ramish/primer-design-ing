from primer_designer.msa_engine import AlignmentMatrix
from primer_designer.conserved_region_detector import find_conserved_regions


def test_find_conserved_region_simple():
    # two identical sequences, length 36 -> should find one conserved region
    seq = 'A'*18 + 'C'*18
    aln = AlignmentMatrix(aligned_sequences=[seq, seq], alignment_length=len(seq), num_sequences=2, full_alignment_text='')
    regions = find_conserved_regions(aln, min_conservation=0.9, max_gap=0.05, min_length=21)
    assert len(regions) >= 1
    r = regions[0]
    assert r.length == 36 or r.length >= 18
    assert r.conservation_score >= 90.0


def test_reject_short_region():
    # conserved region shorter than min_length should be ignored
    seq1 = 'A'*10 + 'C'*10
    seq2 = 'A'*10 + 'C'*10
    aln = AlignmentMatrix(aligned_sequences=[seq1, seq2], alignment_length=len(seq1), num_sequences=2, full_alignment_text='')
    regions = find_conserved_regions(aln, min_conservation=0.9, max_gap=0.05, min_length=21)
    assert len(regions) == 0
