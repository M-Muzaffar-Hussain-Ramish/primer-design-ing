from primer_designer.scoring_engine import score_primer


def test_score_primer_high_quality():
    scored = score_primer(
        primer_id='P-0001',
        conservation_score=95.0,
        tm_basic=62.0,
        tm_advanced=60.0,
        specificity_score=92.0,
        secondary_structure_risk='LOW',
        gc_content=52.0,
        length=20,
    )
    assert scored.composite_score >= 85.0
    assert scored.confidence_rating == 'HIGH'
    assert scored.risk_level == 'LOW'


def test_score_primer_low_quality():
    scored = score_primer(
        primer_id='P-0002',
        conservation_score=60.0,
        tm_basic=45.0,
        tm_advanced=43.0,
        specificity_score=40.0,
        secondary_structure_risk='HIGH',
        gc_content=70.0,
        length=26,
    )
    assert scored.composite_score < 60.0
    assert scored.confidence_rating == 'LOW'
    assert scored.risk_level == 'HIGH'
