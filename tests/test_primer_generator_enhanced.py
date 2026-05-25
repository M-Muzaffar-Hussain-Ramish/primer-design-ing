from primer_designer.primer_generator import generate_candidates_from_region


def test_generate_accepting_candidate():
    # region with balanced composition
    region = ("ATGCGTACGTTAGC" * 6)  # length > 100
    candidates = generate_candidates_from_region(region, length_min=18, length_max=20)
    assert len(candidates) > 0
    # at least one accepted candidate
    assert any(c.is_accepted for c in candidates)


def test_generate_reject_repeat():
    region = ("AT" * 50)  # dinucleotide repeat
    candidates = generate_candidates_from_region(region, length_min=18, length_max=20)
    assert len(candidates) > 0
    # all should be rejected due to dinucleotide repeat
    assert all(not c.is_accepted for c in candidates)
