from dataclasses import dataclass
from typing import List, Tuple, Dict

IUPAC = set("ACGTURYSWKMBDHVNacgturyswkmbdhvn")

@dataclass
class ValidationResult:
    is_valid: bool
    cleaned_sequence: str
    sequence_id: str
    sequence_length: int
    gc_content: float
    composition: Dict[str, int]
    ambiguous_bases: List[Tuple[int, str]]
    errors: List[str]
    warnings: List[str]
    was_rna_converted: bool
    validation_report: str


def _parse_fasta(text: str) -> Tuple[str, str]:
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    if not lines:
        return "", ""
    if lines[0].startswith(">"):
        header = lines[0][1:].strip() or "unknown"
        seq = "".join(lines[1:])
        return header, seq
    # treat as raw
    return "unknown", "".join(lines)


def _compute_gc(seq: str) -> float:
    if not seq:
        return 0.0
    g = seq.count("G") + seq.count("g")
    c = seq.count("C") + seq.count("c")
    return 100.0 * (g + c) / len(seq)


def validate_sequence(raw_input: str, input_format: str = "auto", allow_iupac: bool = True, convert_rna: bool = True) -> ValidationResult:
    errors: List[str] = []
    warnings: List[str] = []
    was_rna_converted = False

    if input_format == "auto":
        seq_id, seq = _parse_fasta(raw_input)
    elif input_format == "fasta":
        seq_id, seq = _parse_fasta(raw_input)
    else:
        seq_id = "unknown"
        seq = raw_input.strip().replace("\n", "").replace(" ", "")

    if convert_rna and ("U" in seq or "u" in seq):
        seq = seq.replace("U", "T").replace("u", "t")
        was_rna_converted = True
        warnings.append("RNA bases converted (U -> T)")

    if len(seq) == 0:
        errors.append("Empty sequence")
        return ValidationResult(False, "", seq_id, 0, 0.0, {}, [], errors, warnings, was_rna_converted, "Empty sequence")

    if len(seq) < 18:
        errors.append("Sequence too short (<18 bp)")

    composition: Dict[str, int] = {}
    ambiguous: List[Tuple[int, str]] = []

    for i, ch in enumerate(seq, start=1):
        composition[ch] = composition.get(ch, 0) + 1
        if ch not in IUPAC:
            errors.append(f"Invalid character at position {i}: {ch}")
        if ch.upper() in "RYSWKMBDHVN":
            ambiguous.append((i, ch))

    gc = _compute_gc(seq)
    if gc < 20 or gc > 80:
        warnings.append(f"Extreme GC content: {gc:.1f}%")

    is_valid = len(errors) == 0

    report_lines = [f"Sequence ID: {seq_id}", f"Length: {len(seq)}", f"GC%: {gc:.1f}"]
    if warnings:
        report_lines.append("Warnings:\n  " + "\n  ".join(warnings))
    if errors:
        report_lines.append("Errors:\n  " + "\n  ".join(errors))

    return ValidationResult(is_valid, seq, seq_id, len(seq), gc, composition, ambiguous, errors, warnings, was_rna_converted, "\n".join(report_lines))
