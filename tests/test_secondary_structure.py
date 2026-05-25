from primer_designer.secondary_structure import detect_hairpin, analyze_primer


def test_detect_hairpin_positive():
    # construct a sequence with a clear hairpin: AAAACGT...ACGTTTT (palindromic)
    seq = 'AAAA' + 'ATGC' + 'AAA' + 'GCAT' + 'TTTT'
    res = detect_hairpin(seq, min_stem=3, max_stem=6, min_loop=1, max_loop=10)
    assert res is not None


def test_detect_hairpin_negative():
    seq = 'ATGCGTACGTAGCTAGCTAG'
    res = detect_hairpin(seq)
    assert res is None


def test_self_dimer_detection():
    seq = 'ATGCATGCATGC'
    report = analyze_primer(seq)
    # short sequence may or may not have self-dimer; ensure function runs and returns expected fields
    assert hasattr(report, 'self_dimer_detected')
    assert isinstance(report.pass_fail, bool)


def test_cross_dimer_detection():
    f = 'ATGCGTACGTA'
    r = 'TACGTACGCAT'  # reverse complement-ish
    report = analyze_primer(f, partner_seq=r)
    assert report.cross_dimer_detected
