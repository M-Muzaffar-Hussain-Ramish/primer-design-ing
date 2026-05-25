import os
from pathlib import Path

import streamlit as st
from primer_designer import PipelineOptions, generate_primer_pair, run_pipeline, validate_sequence


def format_path(output_dir: str, prefix: str, extension: str) -> str:
    return str(Path(output_dir) / f"{prefix}.{extension}")


def display_validation(validation):
    st.markdown("### Validation Summary")
    cols = st.columns(3)
    cols[0].metric("Length", validation.sequence_length)
    cols[1].metric("GC content", f"{validation.gc_content:.1f}%")
    cols[2].metric("Status", "Valid" if validation.is_valid else "Invalid")

    st.markdown("**Sequence ID:** " + validation.sequence_id)
    if validation.was_rna_converted:
        st.warning("RNA bases were converted to DNA (U -> T).")
    if validation.warnings:
        st.warning("Warnings:\n" + "\n".join(validation.warnings))
    if validation.errors:
        st.error("Errors:\n" + "\n".join(validation.errors))

    with st.expander("View full validation report"):
        st.code(validation.validation_report)


def display_pipeline_result(result, output_dir, prefix, include_pdf):
    st.success("Pipeline completed successfully")

    # ── BLAST summary ─────────────────────────────────────────────────────────
    bs = result.get("blast_summary", {})
    if bs:
        cols = st.columns(4)
        cols[0].metric("BLAST hits total", bs.get("total", 0))
        cols[1].metric("Accepted (>= identity threshold)", bs.get("accepted", 0))
        cols[2].metric("Sequences used for MSA", bs.get("sequences_used_for_msa", 0))
        cols[3].metric("Identity threshold", f"{bs.get('identity_threshold', 95):.0f}%")

    if result.get("blast_hits"):
        with st.expander("BLAST hits detail"):
            st.dataframe(result["blast_hits"])

    # ── Conserved region ──────────────────────────────────────────────────────
    with st.expander("Selected conserved region"):
        cr = result.get("selected_region", {})
        rc1, rc2, rc3 = st.columns(3)
        rc1.metric("Length", f"{cr.get('length', 0)} bp")
        rc2.metric("Conservation", f"{cr.get('conservation_score', 0):.1%}")
        rc3.metric("Mean entropy", f"{cr.get('entropy_mean', 0):.3f}")
        st.code(cr.get("consensus_sequence", ""), language=None)

    # ── Primer pair ───────────────────────────────────────────────────────────
    pp = result.get("primer_pair", {})
    fwd = pp.get("forward", {})
    rev = pp.get("reverse", {})

    st.subheader("Designed Primer Pair")

    col_f, col_r = st.columns(2)

    with col_f:
        st.markdown("#### Forward Primer")
        st.code(fwd.get("sequence", "N/A"), language=None)
        m1, m2, m3 = st.columns(3)
        m1.metric("Length", f"{fwd.get('length', 0)} bp")
        m2.metric("Tm", f"{pp.get('tm_forward', 0):.1f} °C")
        m3.metric("GC", f"{fwd.get('gc_content', 0):.1f}%")

    with col_r:
        st.markdown("#### Reverse Primer")
        st.code(rev.get("sequence", "N/A"), language=None)
        m1, m2, m3 = st.columns(3)
        m1.metric("Length", f"{rev.get('length', 0)} bp")
        m2.metric("Tm", f"{pp.get('tm_reverse', 0):.1f} °C")
        m3.metric("GC", f"{rev.get('gc_content', 0):.1f}%")

    comp_cols = st.columns(3)
    tm_diff = pp.get("tm_difference", 0)
    compatible = pp.get("tm_compatible", False)
    comp_cols[0].metric("Tm difference", f"{tm_diff:.2f} °C", delta=None)
    comp_cols[1].metric("Tm compatible", "Yes" if compatible else "No")
    comp_cols[2].metric("Amplicon size", f"{pp.get('amplicon_size_bp', 0)} bp")

    if not compatible:
        st.warning("Tm difference exceeds the configured threshold. Consider relaxing constraints.")

    # ── Exported files ────────────────────────────────────────────────────────
    st.markdown("### Exported Files")
    st.write(
        {
            "JSON": format_path(output_dir, prefix, "json"),
            "TXT": format_path(output_dir, prefix, "txt"),
            "CSV": format_path(output_dir, prefix, "csv"),
            "FASTA": format_path(output_dir, prefix, "fasta"),
            "PDF": format_path(output_dir, prefix, "pdf") if include_pdf else "Not generated",
        }
    )

    with st.expander("MSA alignment output"):
        if result.get("alignment"):
            st.text_area("Full alignment", result["alignment"]["full_alignment_text"], height=240)
        else:
            st.write("No alignment data available.")


def main():
    st.set_page_config(page_title="Primer Designer UI", layout="wide")
    st.title("Primer Designer")
    st.markdown(
        "Design PCR primers, validate nucleotide sequences, and run the full primer-design pipeline with optional export support."
    )

    with st.sidebar:
        st.header("Input")
        input_mode = st.radio("Sequence input mode", ["Raw sequence", "Upload FASTA"])
        sequence_text = ""
        if input_mode == "Raw sequence":
            sequence_text = st.text_area("Paste raw nucleotide sequence", height=240)
        else:
            uploaded = st.file_uploader("Upload a FASTA file", type=["fa", "fasta", "txt"])
            if uploaded is not None:
                sequence_text = uploaded.read().decode("utf-8")

        st.markdown("---")
        st.header("Pipeline settings")
        action = st.selectbox(
            "Action",
            ["Validate sequence", "Generate primers", "Run full pipeline"],
        )

        if action == "Generate primers":
            primer_length = st.number_input(
                "Primer length",
                min_value=18,
                max_value=30,
                value=20,
                help="Length of the forward and reverse primers generated from the sequence.",
            )
        else:
            primer_length = None

        with st.expander("Advanced options", expanded=False):
            use_cache_only = st.checkbox("Use cache-only mode (offline cache only)", value=False)
            disable_ssl_verification = st.checkbox(
                "Disable SSL certificate verification",
                value=False,
                help="Enable this if your network uses a self-signed or corporate TLS proxy certificate.",
            )
            cache_dir = st.text_input("Cache directory", value="cache")
            output_dir = st.text_input("Output directory", value="results")
            prefix = st.text_input("Output filename prefix", value="run")
            include_pdf = st.checkbox("Generate PDF report", value=False)
            primer_length_min = st.number_input(
                "Minimum primer length",
                min_value=16,
                max_value=26,
                value=18,
                help="Minimum primer length for the full pipeline candidate search.",
            )
            primer_length_max = st.number_input(
                "Maximum primer length",
                min_value=18,
                max_value=30,
                value=22,
                help="Maximum primer length for the full pipeline candidate search.",
            )

        run_button = st.button("Execute")

    if not sequence_text.strip():
        st.warning("Enter a sequence or upload a FASTA file to get started.")
        return

    if run_button:
        validation = validate_sequence(sequence_text)
        display_validation(validation)

        if not validation.is_valid:
            return

        if action == "Generate primers":
            if len(validation.cleaned_sequence) < (primer_length or 0) * 2:
                st.error(
                    f"Sequence must be at least {(primer_length or 0) * 2} bp long to generate a forward/reverse pair of {primer_length} bp."
                )
                return

            pair = generate_primer_pair(validation.cleaned_sequence, length=primer_length)
            st.subheader("Primer Pair")
            st.write("**Forward primer**")
            st.code(pair.forward)
            st.write(f"Location: {pair.f_start}-{pair.f_end}")
            st.write("**Reverse primer**")
            st.code(pair.reverse)
            st.write(f"Location: {pair.r_start}-{pair.r_end}")

        elif action == "Run full pipeline":
            opts = PipelineOptions(
                cache_dir=cache_dir,
                output_dir=output_dir,
                prefix=prefix,
                include_pdf=include_pdf,
                use_cache_only=use_cache_only,
                verify_ssl=not disable_ssl_verification,
                primer_length_min=primer_length_min,
                primer_length_max=primer_length_max,
            )
            try:
                with st.spinner("Running the full primer design pipeline..."):
                    result = run_pipeline(validation.cleaned_sequence, options=opts)
                display_pipeline_result(result, output_dir, prefix, include_pdf)
            except Exception as exc:
                st.error(f"Pipeline failed: {exc}")
        else:
            st.info("Validation completed. Choose another action to generate primers or run the pipeline.")


if __name__ == "__main__":
    main()
