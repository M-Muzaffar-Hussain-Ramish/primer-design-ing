"""Primer Designer - minimal viable implementation based on DESIGN.md
"""
from .sequence_validator import ValidationResult, validate_sequence
from .thermodynamics import thermo_profile
from .primer_generator import generate_candidates_from_region, generate_primer_pair, PrimerCandidate, PrimerPair, validate_candidate
from .pipeline_runner import PipelineOptions, run_pipeline

__all__ = [
    "ValidationResult",
    "validate_sequence",
    "thermo_profile",
    "generate_primer_pair",
    "generate_candidates_from_region",
    "validate_candidate",
    "PrimerCandidate",
    "PrimerPair",
    "PipelineOptions",
    "run_pipeline",
]
