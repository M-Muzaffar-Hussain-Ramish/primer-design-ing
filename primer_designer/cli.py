"""Simple CLI for primer_designer"""
import argparse
from primer_designer import validate_sequence, thermo_profile, generate_primer_pair


def main():
    p = argparse.ArgumentParser(prog="primer_designer")
    p.add_argument("input", help="FASTA file or raw sequence string")
    p.add_argument("--length", type=int, default=20, help="Primer length")
    args = p.parse_args()

    # Read input (file if path exists)
    try:
        with open(args.input, "r", encoding="utf-8") as fh:
            raw = fh.read()
    except FileNotFoundError:
        raw = args.input

    v = validate_sequence(raw)
    if not v.is_valid:
        print("Validation failed:\n", v.validation_report)
        return

    print(f"Sequence validated: {v.sequence_id} | Length: {v.sequence_length} | GC%: {v.gc_content:.1f}%")
    pair = generate_primer_pair(v.cleaned_sequence, length=args.length)
    tf = thermo_profile(pair.forward)
    tr = thermo_profile(pair.reverse)

    print("\n=== Selected Primers ===")
    print(f"Forward ({pair.f_start}-{pair.f_end}): {pair.forward} | Tm_basic={tf.tm_basic:.1f} | GC={tf.gc_percent:.1f}%")
    print(f"Reverse ({pair.r_start}-{pair.r_end}): {pair.reverse} | Tm_basic={tr.tm_basic:.1f} | GC={tr.gc_percent:.1f}%")


if __name__ == "__main__":
    main()
