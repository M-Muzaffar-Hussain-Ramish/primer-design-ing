from primer_designer.primer_generator import generate_primer_pair


def test_generate_pair():
    seq = "A" * 50 + "C" * 50
    pair = generate_primer_pair(seq, length=20)
    assert pair.forward == seq[0:20]
    # reverse primer should be reverse-complement of last 20
    last20 = seq[-20:]
    # complement of C is G
    assert pair.reverse == last20[::-1].translate(str.maketrans('ACGT', 'TGCA'))
