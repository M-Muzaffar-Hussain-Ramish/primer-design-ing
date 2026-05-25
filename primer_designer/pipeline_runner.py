import json
import os
from dataclasses import asdict, dataclass
from typing import Any, Dict, List, Optional

from primer_designer.blast_engine import BlastConfig, BlastHit, run_blast
from primer_designer.conserved_region_detector import find_conserved_regions
from primer_designer.export_manager import ExportOptions, export_csv, export_fasta, export_json, export_pdf, export_txt
from primer_designer.msa_engine import MSAConfig, align_sequences
from primer_designer.primer_generator import generate_candidates_from_region
from primer_designer.scoring_engine import score_primer
from primer_designer.secondary_structure import analyze_primer
from primer_designer.sequence_validator import validate_sequence
from primer_designer.specificity_checker import run_primer_specificity


@dataclass
class PipelineOptions:
    cache_dir: str = "cache"
    output_dir: str = "results"
    prefix: str = "pipeline_report"
    include_pdf: bool = False
    primer_length_min: int = 18
    primer_length_max: int = 22
    min_conserved_length: int = 18
    use_cache_only: bool = True


def _serialize_hits(hits: List[BlastHit]) -> List[Dict[str, Any]]:
    return [asdict(hit) for hit in hits]


def run_pipeline(raw_input: str, options: Optional[PipelineOptions] = None) -> Dict[str, Any]:
    if options is None:
        options = PipelineOptions()

    validation = validate_sequence(raw_input)
    if not validation.is_valid:
        raise ValueError("Invalid input sequence: " + "; ".join(validation.errors))

    blast_cfg = BlastConfig(use_cache=True, cache_only=options.use_cache_only)
    blast_hits = run_blast(validation.cleaned_sequence, config=blast_cfg, cache_dir=options.cache_dir)

    msa_cfg = MSAConfig(use_cache=options.use_cache_only)
    alignment = align_sequences([validation.cleaned_sequence, validation.cleaned_sequence], config=msa_cfg, cache_dir=options.cache_dir)

    regions = find_conserved_regions(alignment, min_length=options.min_conserved_length)
    if not regions:
        raise ValueError("No conserved regions detected")
    region = regions[0]

    candidates = generate_candidates_from_region(region.consensus_sequence, length_min=options.primer_length_min, length_max=options.primer_length_max)
    if not candidates:
        raise ValueError("No primer candidates generated")

    scored_candidates = []
    for candidate in candidates:
        structure = analyze_primer(candidate.sequence)
        try:
            specificity = run_primer_specificity(candidate.sequence, candidate.candidate_id, cache_dir=options.cache_dir)
            spec_score = specificity.specificity_score
            risk = specificity.risk_level
        except RuntimeError:
            spec_score = 100.0
            risk = "LOW"
        score = score_primer(
            primer_id=candidate.candidate_id,
            conservation_score=region.conservation_score,
            tm_basic=candidate.tm_basic,
            tm_advanced=candidate.tm_advanced,
            specificity_score=spec_score,
            secondary_structure_risk=structure.overall_risk,
            gc_content=candidate.gc_content,
            length=candidate.length,
        )
        scored_candidates.append({
            **asdict(candidate),
            "structure": asdict(structure),
            "specificity": {
                "score": spec_score,
                "risk_level": risk,
            },
            "scoring": asdict(score),
        })

    selected = sorted(scored_candidates, key=lambda x: x["scoring"]["composite_score"], reverse=True)[0]

    primer_table = [
        {
            "id": c["candidate_id"],
            "sequence": c["sequence"],
            "composite_score": c["scoring"]["composite_score"],
            "risk_level": c["scoring"]["risk_level"],
        }
        for c in scored_candidates
    ]

    result = {
        "validation": asdict(validation),
        "blast_hits": _serialize_hits(blast_hits),
        "alignment": {
            "alignment_length": alignment.alignment_length,
            "num_sequences": alignment.num_sequences,
            "full_alignment_text": alignment.full_alignment_text,
        },
        "selected_region": {
            "region_id": region.region_id,
            "start_position": region.start_position,
            "end_position": region.end_position,
            "length": region.length,
            "conservation_score": region.conservation_score,
            "consensus_sequence": region.consensus_sequence,
        },
        "candidate_count": len(scored_candidates),
        "selected_candidate": selected,
    }

    export_json(result, ExportOptions(output_dir=options.output_dir, prefix=options.prefix))
    export_txt(result, ExportOptions(output_dir=options.output_dir, prefix=options.prefix))
    export_csv(primer_table, ExportOptions(output_dir=options.output_dir, prefix=options.prefix))
    export_fasta(
        [{"id": c["candidate_id"], "sequence": c["sequence"]} for c in scored_candidates],
        ExportOptions(output_dir=options.output_dir, prefix=options.prefix),
    )
    export_pdf(result, ExportOptions(output_dir=options.output_dir, prefix=options.prefix, include_pdf=options.include_pdf))

    return result
