from primer_designer.thermodynamics import thermo_profile


def test_thermo_basic_and_advanced():
    seq = "ATGCTAGCTAGCTAGCTAGC"
    t = thermo_profile(seq)
    assert t.length == len(seq)
    assert abs(t.tm_basic - (2*(seq.count('A')+seq.count('T')) + 4*(seq.count('G')+seq.count('C')))) < 1e-6
    assert isinstance(t.tm_advanced, float)
    assert 0 <= t.gc_percent <= 100
