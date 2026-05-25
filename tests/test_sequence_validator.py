import pytest
from primer_designer.sequence_validator import validate_sequence


def test_validate_raw_sequence():
    seq = "ATGCGTACGTTAGCCTAGCT"
    r = validate_sequence(seq, input_format="raw")
    assert r.is_valid
    assert r.sequence_length == len(seq)
    assert not r.was_rna_converted


def test_rna_conversion():
    seq = "AUGCGUACGUU"
    r = validate_sequence(seq, input_format="raw")
    assert r.was_rna_converted
    assert "U" not in r.cleaned_sequence.upper()


def test_invalid_character():
    seq = "ATGXB"
    r = validate_sequence(seq, input_format="raw")
    assert not r.is_valid
    assert any("Invalid character" in e for e in r.errors)
