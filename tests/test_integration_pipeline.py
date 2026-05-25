import json
from pathlib import Path

from primer_designer.blast_engine import BlastConfig, _query_hash
from primer_designer.msa_engine import MSAConfig, _sequences_hash
from primer_designer.pipeline_runner import PipelineOptions, run_pipeline


def write_blast_cache(cache_dir: Path, sequence: str):
    cfg = BlastConfig(program="blastn", database="nt", identity_threshold=95.0, coverage_threshold=70.0, max_hits=20, e_value_cutoff=1e-10, use_cache=True, cache_only=True)
    qh = _query_hash(sequence, cfg)
    path = cache_dir / f"{qh}.json"
    sample = {
        "query_hash": qh,
        "hits": [
            {
                "accession": "NM_000000",
                "organism": "Homo sapiens",
                "alignment_length": len(sequence),
                "query_coverage": 100.0,
                "identity_percent": 100.0,
                "mismatches": 0,
                "gaps": 0,
                "e_value": 0.0,
                "bit_score": 200.0,
                "aligned_segment_preview": sequence[:80],
                "query_start": 1,
                "query_end": len(sequence),
                "subject_start": 1,
                "subject_end": len(sequence),
                "strand": "plus/plus",
                "is_accepted": True,
                "rejection_reason": None,
                "raw_score": 200.0,
                "subject_sequence": sequence,
            }
        ],
    }
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(sample, fh)


def write_msa_cache(cache_dir: Path, sequence: str):
    cfg = MSAConfig(use_cache=True)
    qh = _sequences_hash([sequence, sequence], cfg)
    path = cache_dir / f"msa_{qh}.json"
    sample = {
        "aligned_sequences": [sequence, sequence],
        "alignment_length": len(sequence),
        "num_sequences": 2,
        "full_alignment_text": f"CLUSTAL W\nseq_1    {sequence}\nseq_2    {sequence}\n",
    }
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(sample, fh)


def test_full_pipeline_with_caches(tmp_path):
    raw = "ATGCGTACGTTAGCCTAGCT"
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    write_blast_cache(cache_dir, raw)
    write_msa_cache(cache_dir, raw)

    output_dir = tmp_path / "results"
    options = PipelineOptions(cache_dir=str(cache_dir), output_dir=str(output_dir), prefix="integration", include_pdf=False, use_cache_only=True)
    result = run_pipeline(raw, options=options)

    assert result["validation"]["is_valid"]
    assert result["blast_hits"]
    assert result["alignment"]["alignment_length"] == len(raw)
    assert result["selected_region"]["length"] >= 18
    assert "selected_candidate" in result
    assert Path(output_dir / "integration.json").exists()
    assert Path(output_dir / "integration.txt").exists()
    assert Path(output_dir / "integration.csv").exists()
    assert Path(output_dir / "integration.fasta").exists()
