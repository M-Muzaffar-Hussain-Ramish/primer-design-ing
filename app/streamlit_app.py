import streamlit as st
from primer_designer import validate_sequence, generate_primer_pair, run_pipeline, PipelineOptions


def main():
    st.set_page_config(page_title="Primer Designer UI", layout="wide")
    st.title("Primer Designer")
    st.markdown(
        "Use this interface to validate sequences, generate primer pairs, or run the full pipeline using cached BLAST/MSA results."
    )

    with st.sidebar:
        st.header("Input")
        input_mode = st.radio("Sequence input mode", ["Raw sequence", "Upload FASTA"])
        sequence_text = ""
        if input_mode == "Raw sequence":
            sequence_text = st.text_area("Paste raw nucleotide sequence", height=200)
        else:
            uploaded = st.file_uploader("Upload a FASTA file", type=["fa", "fasta", "txt"])
            if uploaded is not None:
                sequence_text = uploaded.read().decode("utf-8")
        st.markdown("---")
        st.header("Parameters")
        primer_length = st.number_input("Primer length", min_value=18, max_value=25, value=20)
        use_cache_only = st.checkbox("Use cache-only mode", value=True)
        cache_dir = st.text_input("Cache directory", value="cache")
        output_dir = st.text_input("Output directory", value="results")
        prefix = st.text_input("Output filename prefix", value="run")
        include_pdf = st.checkbox("Generate PDF report", value=False)
        action = st.selectbox(
            "Action",
            ["Validate sequence", "Generate primers", "Run full pipeline"],
        )
        run_button = st.button("Execute")

    if not sequence_text:
        st.warning("Enter a sequence or upload a FASTA file to get started.")
        return

    if run_button:
        validation = validate_sequence(sequence_text)
        st.subheader("Validation Result")
        if validation.is_valid:
            st.success("Sequence is valid")
        else:
            st.error("Sequence is invalid")
        st.text(validation.validation_report)

        if not validation.is_valid:
            return

        if action == "Generate primers":
            if len(validation.cleaned_sequence) < primer_length * 2:
                st.error(
                    f"Sequence must be at least {primer_length * 2} bp long to generate a forward/reverse pair of {primer_length} bp."
                )
                return
            pair = generate_primer_pair(validation.cleaned_sequence, length=primer_length)
            st.subheader("Primer Pair")
            st.write("Forward primer:")
            st.code(pair.forward)
            st.write(f"Location: {pair.f_start}-{pair.f_end}")
            st.write("Reverse primer:")
            st.code(pair.reverse)
            st.write(f"Location: {pair.r_start}-{pair.r_end}")

        elif action == "Run full pipeline":
            st.subheader("Pipeline Execution")
            opts = PipelineOptions(
                cache_dir=cache_dir,
                output_dir=output_dir,
                prefix=prefix,
                include_pdf=include_pdf,
                use_cache_only=use_cache_only,
            )
            try:
                result = run_pipeline(sequence_text, options=opts)
                st.success("Pipeline completed successfully")
                st.subheader("Selected Candidate")
                st.json(result["selected_candidate"])
                st.subheader("Exported Files")
                st.write(
                    {
                        "JSON": f"{output_dir}/{prefix}.json",
                        "TXT": f"{output_dir}/{prefix}.txt",
                        "CSV": f"{output_dir}/{prefix}.csv",
                        "FASTA": f"{output_dir}/{prefix}.fasta",
                        "PDF": f"{output_dir}/{prefix}.pdf" if include_pdf else "Not generated",
                    }
                )
                st.subheader("Pipeline Result Summary")
                st.json({
                    "selected_region": result["selected_region"],
                    "candidate_count": result["candidate_count"],
                })
            except Exception as exc:
                st.error(f"Pipeline failed: {exc}")

        else:
            st.info("Validation completed. Select another action to continue.")


if __name__ == "__main__":
    main()
