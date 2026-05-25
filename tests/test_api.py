"""Tests for the FastAPI REST API endpoints."""
import json
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from api.main import app
from primer_designer.blast_engine import BlastConfig, _query_hash
from primer_designer.msa_engine import MSAConfig, _sequences_hash

client = TestClient(app)

SEQ = "ATGCGTACGTTAGCCTAGCTATGCGTACGTTAGCCTAGCT"  # 40 bp, valid DNA


# ---------------------------------------------------------------------------
# /health
# ---------------------------------------------------------------------------

def test_health():
    resp = client.get("/health")
    assert resp.status_code == 200
    assert resp.json()["status"] == "ok"


# ---------------------------------------------------------------------------
# /api/v1/validate
# ---------------------------------------------------------------------------

def test_validate_valid_sequence():
    resp = client.post("/api/v1/validate", json={"sequence": SEQ})
    assert resp.status_code == 200
    data = resp.json()
    assert data["is_valid"] is True
    assert data["sequence_length"] == len(SEQ)
    assert 0 < data["gc_content"] < 100


def test_validate_rna_conversion():
    rna = "AUGCGUACGU" * 3  # 30 bases RNA
    resp = client.post("/api/v1/validate", json={"sequence": rna, "convert_rna": True})
    assert resp.status_code == 200
    data = resp.json()
    assert data["is_valid"] is True
    assert data["was_rna_converted"] is True


def test_validate_invalid_chars():
    resp = client.post("/api/v1/validate", json={"sequence": "ATGCXYZ123"})
    assert resp.status_code == 200
    data = resp.json()
    assert data["is_valid"] is False
    assert data["errors"]


def test_validate_too_short():
    resp = client.post("/api/v1/validate", json={"sequence": "ATGC"})
    assert resp.status_code == 200
    data = resp.json()
    assert data["is_valid"] is False


# ---------------------------------------------------------------------------
# /api/v1/thermodynamics
# ---------------------------------------------------------------------------

def test_thermodynamics():
    primer = "ATGCGTACGTTAGCCTAGCT"  # 20 bp
    resp = client.post("/api/v1/thermodynamics", json={"sequence": primer})
    assert resp.status_code == 200
    data = resp.json()
    assert data["length"] == 20
    assert 0 < data["gc_percent"] <= 100
    assert data["tm_basic_wallace"] > 0
    assert data["tm_advanced"] > 0
    assert abs(data["tm_mean"] - (data["tm_basic_wallace"] + data["tm_advanced"]) / 2) < 0.01


# ---------------------------------------------------------------------------
# /api/v1/secondary-structure
# ---------------------------------------------------------------------------

def test_secondary_structure_no_hairpin():
    # A random sequence unlikely to form strong structures
    resp = client.post("/api/v1/secondary-structure", json={"sequence": "ATGCATGCATGCATGCATGC"})
    assert resp.status_code == 200
    data = resp.json()
    assert "overall_risk" in data
    assert data["overall_risk"] in ("LOW", "MEDIUM", "HIGH", "FAIL")


def test_secondary_structure_with_partner():
    resp = client.post(
        "/api/v1/secondary-structure",
        json={"sequence": "ATGCATGCATGCATGCATGC", "partner_sequence": "GCATGCATGCATGCATGCAT"},
    )
    assert resp.status_code == 200
    data = resp.json()
    assert "cross_dimer_detected" in data


# ---------------------------------------------------------------------------
# /api/v1/primer-pair
# ---------------------------------------------------------------------------

def test_primer_pair():
    resp = client.post("/api/v1/primer-pair", json={"sequence": SEQ, "length": 18})
    assert resp.status_code == 200
    data = resp.json()
    assert len(data["forward"]["sequence"]) == 18
    assert len(data["reverse"]["sequence"]) == 18
    assert data["amplicon_size"] > 0
    assert "tm_compatible" in data


def test_primer_pair_too_short_sequence():
    resp = client.post("/api/v1/primer-pair", json={"sequence": "ATGCATGC" * 4, "length": 20})
    # 32 bp < 40 bp required → 422
    assert resp.status_code == 422


# ---------------------------------------------------------------------------
# /api/v1/candidates
# ---------------------------------------------------------------------------

def test_candidates():
    region = "ATGCGTACGTTAGCCTAGCTATGCGTACGTTAGCCTAGCT"
    resp = client.post(
        "/api/v1/candidates",
        json={"region_sequence": region, "length_min": 18, "length_max": 20},
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["total_candidates"] > 0
    assert isinstance(data["candidates"], list)
    c = data["candidates"][0]
    assert "sequence" in c
    assert "gc_content" in c
    assert "is_accepted" in c


# ---------------------------------------------------------------------------
# /api/v1/pipeline  (cache-only mode using fixture caches)
# ---------------------------------------------------------------------------

def _write_blast_cache(cache_dir: Path, sequence: str):
    cfg = BlastConfig(use_cache=True, cache_only=True, identity_threshold=95.0, coverage_threshold=70.0, max_hits=20)
    qh = _query_hash(sequence, cfg)
    path = cache_dir / f"{qh}.json"
    payload = {
        "query_hash": qh,
        "hits": [
            {
                "accession": "NM_000001",
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
    path.write_text(json.dumps(payload))


def _write_msa_cache(cache_dir: Path, sequences: list):
    cfg = MSAConfig(use_cache=True)
    qh = _sequences_hash(sequences, cfg)
    path = cache_dir / f"msa_{qh}.json"
    aln = sequences[0]
    payload = {
        "aligned_sequences": sequences,
        "alignment_length": len(aln),
        "num_sequences": len(sequences),
        "full_alignment_text": "CLUSTAL W\n" + "\n".join(
            f"seq_{i+1}    {s}" for i, s in enumerate(sequences)
        ),
    }
    path.write_text(json.dumps(payload))


def test_pipeline_cache_only(tmp_path):
    seq = "ATGCGTACGTTAGCCTAGCTATGCGTACGTTAGCCTAGCT"  # 40 bp
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    results_dir = tmp_path / "results"

    _write_blast_cache(cache_dir, seq)
    _write_msa_cache(cache_dir, [seq, seq])

    resp = client.post(
        "/api/v1/pipeline",
        json={
            "sequence": seq,
            "cache_dir": str(cache_dir),
            "output_dir": str(results_dir),
            "prefix": "api_test",
            "use_cache_only": True,
            "primer_length_min": 18,
            "primer_length_max": 20,
        },
    )
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["validation"]["is_valid"]
    assert "selected_candidate" in data
    assert (results_dir / "api_test.json").exists()


def test_pipeline_invalid_sequence():
    resp = client.post(
        "/api/v1/pipeline",
        json={"sequence": "ATGCXYZ!", "use_cache_only": True},
    )
    assert resp.status_code == 422
