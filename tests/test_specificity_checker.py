import json
from pathlib import Path
from primer_designer.specificity_checker import run_primer_specificity
from primer_designer.blast_engine import _query_hash, BlastConfig


def test_specificity_report_from_cache(tmp_path):
    primer = "ATGCGTACGTTAGCCTAGCT"
    cfg = BlastConfig(program="blastn", database="nt", identity_threshold=70.0, coverage_threshold=0.0, max_hits=20, e_value_cutoff=1e-3, use_cache=True)
    qh = _query_hash(primer, cfg)
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    path = cache_dir / f"{qh}.json"
    sample = {
        "query_hash": qh,
        "hits": [
            {
                "accession": "NM_000000",
                "organism": "Homo sapiens",
                "alignment_length": 20,
                "query_coverage": 100.0,
                "identity_percent": 100.0,
                "mismatches": 0,
                "gaps": 0,
                "e_value": 0.0,
                "bit_score": 200.0,
                "aligned_segment_preview": "ATGCGT...",
                "query_start": 1,
                "query_end": 20,
                "subject_start": 100,
                "subject_end": 119,
                "strand": "plus/plus",
                "is_accepted": True,
                "rejection_reason": None,
                "raw_score": 200.0,
            }
        ],
    }
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(sample, fh)

    report = run_primer_specificity(primer, "PRIMER-1", cache_dir=str(cache_dir))
    assert report.primer_id == "PRIMER-1"
    assert report.total_hits == 1
    assert report.risk_level == "LOW"
    assert report.recommendation == "ACCEPTED"


def test_specificity_report_pseudogene_risk(tmp_path):
    primer = "ATGCGTACGTTAGCCTAGCT"
    cfg = BlastConfig(program="blastn", database="nt", identity_threshold=70.0, coverage_threshold=0.0, max_hits=20, e_value_cutoff=1e-3, use_cache=True)
    qh = _query_hash(primer, cfg)
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    path = cache_dir / f"{qh}.json"
    sample = {
        "query_hash": qh,
        "hits": [
            {
                "accession": "XYZ_000001",
                "organism": "Human pseudogene",
                "alignment_length": 20,
                "query_coverage": 100.0,
                "identity_percent": 96.0,
                "mismatches": 0,
                "gaps": 0,
                "e_value": 0.0,
                "bit_score": 190.0,
                "aligned_segment_preview": "ATGCGT...",
                "query_start": 1,
                "query_end": 20,
                "subject_start": 101,
                "subject_end": 120,
                "strand": "plus/plus",
                "is_accepted": True,
                "rejection_reason": None,
                "raw_score": 190.0,
            }
        ],
    }
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(sample, fh)

    report = run_primer_specificity(primer, "PRIMER-2", cache_dir=str(cache_dir))
    assert report.risk_level == "CRITICAL"
    assert report.recommendation == "REJECTED"
