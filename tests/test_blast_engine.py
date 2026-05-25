import os
import json
from primer_designer.blast_engine import run_blast, BlastConfig, _query_hash, _cache_path


def test_load_from_cache(tmp_path, monkeypatch):
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    query = "ATGCGTACGTTAGCCTAGCT"
    cfg = BlastConfig(use_cache=True)
    qh = _query_hash(query, cfg)
    path = cache_dir / f"{qh}.json"
    sample = {"query_hash": qh, "hits": [
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
            "raw_score": 200.0
        }
    ]}
    with open(path, 'w', encoding='utf-8') as fh:
        json.dump(sample, fh)

    hits = run_blast(query, config=cfg, cache_dir=str(cache_dir))
    assert len(hits) == 1
    assert hits[0].accession == "NM_000000"
    assert hits[0].organism.startswith("Homo sapiens")


def test_no_cache_and_no_biopython(tmp_path, monkeypatch):
    # if Biopython is not available, run_blast should raise informative error
    cache_dir = tmp_path / "cache2"
    cache_dir.mkdir()
    query = "ATGC"
    cfg = BlastConfig(use_cache=False)

    # ensure Biopython import fails
    import sys
    monkeypatch.setitem(sys.modules, 'Bio', None)
    # Remove Bio if present
    sys.modules.pop('Bio', None)

    try:
        run_blast(query, config=cfg, cache_dir=str(cache_dir))
        raised = False
    except RuntimeError as e:
        raised = True
        assert 'Biopython' in str(e) or 'cached' in str(e)
    assert raised
