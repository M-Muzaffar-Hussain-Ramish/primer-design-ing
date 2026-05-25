"""Primer Designer REST API — powered by FastAPI.

Start with:
    uvicorn api.main:app --host 0.0.0.0 --port 8000 --reload

Set environment variables (see .env.example) before starting.
"""
from __future__ import annotations

import os
from typing import Any, Dict, List, Optional

from fastapi import BackgroundTasks, FastAPI, HTTPException, UploadFile, File
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field, field_validator

from primer_designer.blast_engine import BlastConfig, run_blast
from primer_designer.pipeline_runner import PipelineOptions, run_pipeline
from primer_designer.primer_generator import generate_candidates_from_region, generate_primer_pair
from primer_designer.scoring_engine import score_primer
from primer_designer.secondary_structure import analyze_primer
from primer_designer.sequence_validator import validate_sequence
from primer_designer.thermodynamics import thermo_profile

# ---------------------------------------------------------------------------
# App setup
# ---------------------------------------------------------------------------

app = FastAPI(
    title="Primer Designer API",
    description=(
        "Research-grade PCR primer design pipeline. "
        "Validates sequences, runs BLAST homology searches, "
        "performs global MSA, detects conserved regions, "
        "generates and scores primer candidates."
    ),
    version="1.0.0",
    docs_url="/docs",
    redoc_url="/redoc",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=os.environ.get("CORS_ORIGINS", "*").split(","),
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

CACHE_DIR = os.environ.get("CACHE_DIR", "cache")
OUTPUT_DIR = os.environ.get("OUTPUT_DIR", "results")


# ---------------------------------------------------------------------------
# Request / Response models
# ---------------------------------------------------------------------------

class SequenceRequest(BaseModel):
    sequence: str = Field(..., min_length=1, description="Raw DNA/RNA sequence or FASTA text")
    allow_iupac: bool = Field(True, description="Allow IUPAC ambiguity codes")
    convert_rna: bool = Field(True, description="Convert U → T if RNA detected")


class PrimerPairRequest(BaseModel):
    sequence: str = Field(..., min_length=36, description="DNA sequence (≥2× primer_length)")
    length: int = Field(20, ge=15, le=35, description="Primer length in bp")
    offset_from_ends: int = Field(0, ge=0, description="Bases to skip from each end")


class ThermodynamicsRequest(BaseModel):
    sequence: str = Field(..., min_length=10, description="Primer sequence")


class StructureRequest(BaseModel):
    sequence: str = Field(..., min_length=10, description="Primer sequence")
    partner_sequence: Optional[str] = Field(None, description="Partner primer for cross-dimer check")


class CandidatesRequest(BaseModel):
    region_sequence: str = Field(..., min_length=18, description="Conserved region consensus sequence")
    length_min: int = Field(18, ge=15, le=30)
    length_max: int = Field(22, ge=15, le=35)


class PipelineRequest(BaseModel):
    sequence: str = Field(..., min_length=18, description="Input DNA sequence or FASTA")
    cache_dir: str = Field("cache", description="Cache directory for BLAST/MSA results")
    output_dir: str = Field("results", description="Directory for exported result files")
    prefix: str = Field("run", description="Filename prefix for exported files")
    include_pdf: bool = Field(False, description="Generate PDF report")
    primer_length_min: int = Field(18, ge=15, le=30)
    primer_length_max: int = Field(22, ge=15, le=35)
    min_conserved_length: int = Field(18, ge=10, le=100)
    use_cache_only: bool = Field(
        False,
        description="Only use pre-computed caches — no network BLAST calls",
    )
    verify_ssl: bool = Field(True, description="Verify HTTPS certificates for NCBI requests")

    @field_validator("primer_length_max")
    @classmethod
    def max_gte_min(cls, v: int, info: Any) -> int:
        min_val = info.data.get("primer_length_min", 18)
        if v < min_val:
            raise ValueError("primer_length_max must be >= primer_length_min")
        return v


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------

@app.get("/health", tags=["Health"])
def health_check() -> Dict[str, str]:
    """Liveness probe — always returns 200 OK when the service is up."""
    return {"status": "ok", "version": app.version}


@app.post("/api/v1/validate", tags=["Sequence"])
def validate(req: SequenceRequest) -> Dict[str, Any]:
    """Validate a nucleotide sequence and return GC%, composition, and any errors."""
    result = validate_sequence(
        req.sequence,
        allow_iupac=req.allow_iupac,
        convert_rna=req.convert_rna,
    )
    return {
        "is_valid": result.is_valid,
        "sequence_id": result.sequence_id,
        "sequence_length": result.sequence_length,
        "gc_content": result.gc_content,
        "composition": result.composition,
        "ambiguous_bases": result.ambiguous_bases,
        "was_rna_converted": result.was_rna_converted,
        "errors": result.errors,
        "warnings": result.warnings,
        "validation_report": result.validation_report,
    }


@app.post("/api/v1/validate/fasta", tags=["Sequence"])
async def validate_fasta(file: UploadFile = File(...)) -> Dict[str, Any]:
    """Upload a FASTA file and validate the first sequence."""
    contents = await file.read()
    raw = contents.decode("utf-8", errors="replace")
    result = validate_sequence(raw)
    if not result.is_valid:
        raise HTTPException(status_code=422, detail=result.errors)
    return {
        "is_valid": result.is_valid,
        "sequence_id": result.sequence_id,
        "sequence_length": result.sequence_length,
        "gc_content": result.gc_content,
        "warnings": result.warnings,
    }


@app.post("/api/v1/thermodynamics", tags=["Primer"])
def thermodynamics(req: ThermodynamicsRequest) -> Dict[str, Any]:
    """Calculate Tm (Wallace rule & advanced formula) and GC% for a primer sequence."""
    seq = req.sequence.strip().upper()
    tp = thermo_profile(seq)
    return {
        "sequence": seq,
        "length": tp.length,
        "gc_percent": tp.gc_percent,
        "tm_basic_wallace": tp.tm_basic,
        "tm_advanced": tp.tm_advanced,
        "tm_mean": (tp.tm_basic + tp.tm_advanced) / 2.0,
        "count_A": tp.count_A,
        "count_T": tp.count_T,
        "count_G": tp.count_G,
        "count_C": tp.count_C,
    }


@app.post("/api/v1/secondary-structure", tags=["Primer"])
def secondary_structure(req: StructureRequest) -> Dict[str, Any]:
    """Analyse a primer for hairpin, self-dimer, and cross-dimer secondary structures."""
    report = analyze_primer(req.sequence, req.partner_sequence)
    return {
        "primer_sequence": req.sequence,
        "partner_sequence": req.partner_sequence,
        "hairpin_detected": report.hairpin_detected,
        "hairpin_delta_g": report.hairpin_delta_g,
        "self_dimer_detected": report.self_dimer_detected,
        "self_dimer_delta_g": report.self_dimer_delta_g,
        "cross_dimer_detected": report.cross_dimer_detected,
        "cross_dimer_delta_g": report.cross_dimer_delta_g,
        "three_prime_binding": report.three_prime_binding,
        "overall_risk": report.overall_risk,
        "pass_fail": report.pass_fail,
        "failure_explanation": report.failure_explanation,
    }


@app.post("/api/v1/primer-pair", tags=["Primer"])
def primer_pair(req: PrimerPairRequest) -> Dict[str, Any]:
    """Generate a simple forward/reverse primer pair from the 5' and 3' ends of the sequence."""
    validation = validate_sequence(req.sequence)
    if not validation.is_valid:
        raise HTTPException(status_code=422, detail=validation.errors)

    seq = validation.cleaned_sequence
    if len(seq) < req.length * 2 + req.offset_from_ends:
        raise HTTPException(
            status_code=422,
            detail=f"Sequence too short. Need ≥{req.length * 2 + req.offset_from_ends} bp, got {len(seq)} bp.",
        )

    pair = generate_primer_pair(seq, length=req.length, offset_from_ends=req.offset_from_ends)
    tf = thermo_profile(pair.forward)
    tr = thermo_profile(pair.reverse)
    tm_diff = abs(tf.tm_advanced - tr.tm_advanced)

    return {
        "forward": {
            "sequence": pair.forward,
            "start": pair.f_start,
            "end": pair.f_end,
            "tm_basic": tf.tm_basic,
            "tm_advanced": tf.tm_advanced,
            "gc_percent": tf.gc_percent,
        },
        "reverse": {
            "sequence": pair.reverse,
            "start": pair.r_start,
            "end": pair.r_end,
            "tm_basic": tr.tm_basic,
            "tm_advanced": tr.tm_advanced,
            "gc_percent": tr.gc_percent,
        },
        "amplicon_size": pair.r_end - pair.f_start + 1,
        "tm_difference": round(tm_diff, 2),
        "tm_compatible": tm_diff <= 2.0,
    }


@app.post("/api/v1/candidates", tags=["Primer"])
def candidates(req: CandidatesRequest) -> Dict[str, Any]:
    """Enumerate all primer candidates from a conserved region sequence."""
    all_candidates = generate_candidates_from_region(
        req.region_sequence,
        length_min=req.length_min,
        length_max=req.length_max,
    )
    accepted = [c for c in all_candidates if c.is_accepted]
    return {
        "total_candidates": len(all_candidates),
        "accepted_candidates": len(accepted),
        "candidates": [
            {
                "id": c.candidate_id,
                "sequence": c.sequence,
                "length": c.length,
                "gc_content": c.gc_content,
                "tm_basic": c.tm_basic,
                "tm_advanced": c.tm_advanced,
                "gc_clamp": c.gc_clamp,
                "is_accepted": c.is_accepted,
                "rejection_reasons": c.rejection_reasons,
            }
            for c in all_candidates
        ],
    }


@app.post("/api/v1/pipeline", tags=["Pipeline"])
def pipeline(req: PipelineRequest) -> Dict[str, Any]:
    """Run the full primer design pipeline.

    Steps: validate → BLAST → MSA → conserved regions → candidates → score → export.
    Use ``use_cache_only=true`` to run entirely offline with pre-populated caches.

    Returns the full result dictionary including the best-scoring primer candidate
    and paths to all exported files.
    """
    validation = validate_sequence(req.sequence)
    if not validation.is_valid:
        raise HTTPException(status_code=422, detail=validation.errors)

    opts = PipelineOptions(
        cache_dir=req.cache_dir,
        output_dir=req.output_dir,
        prefix=req.prefix,
        include_pdf=req.include_pdf,
        primer_length_min=req.primer_length_min,
        primer_length_max=req.primer_length_max,
        min_conserved_length=req.min_conserved_length,
        use_cache_only=req.use_cache_only,
        verify_ssl=req.verify_ssl,
    )

    try:
        result = run_pipeline(validation.cleaned_sequence, options=opts)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))
    except RuntimeError as exc:
        raise HTTPException(status_code=503, detail=str(exc))

    return result
