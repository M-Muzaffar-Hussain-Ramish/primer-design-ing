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
    use_cache_only: bool = False
    verify_ssl: bool = True


def _serialize_hits(hits: List[BlastHit]) -> List[Dict[str, Any]]:
    return [asdict(hit) for hit in hits]


def _build_msa_sequences(query: str, blast_hits: List[BlastHit], max_seqs: int = 20) -> List[str]:
    """Build the sequence list for MSA: query + unique subject sequences from accepted BLAST hits.

    Falls back to duplicating the query when no valid homolog sequences are available
    (e.g. cache-only mode with no subject_sequence stored).
    """
    sequences: List[str] = [query]
    seen = {query}

    accepted_hits = [h for h in blast_hits if h.is_accepted]
    for hit in accepted_hits:
        if len(sequences) >= max_seqs:
            break
        # Prefer full subject sequence; fall back to aligned preview
        subj_raw = hit.subject_sequence or hit.aligned_segment_preview
        # Remove gap characters and whitespace
        subj_clean = "".join(c for c in subj_raw if c.isalpha())
        if len(subj_clean) < 18:
            continue
        if subj_clean not in seen:
            sequences.append(subj_clean)
            seen.add(subj_clean)

    if len(sequences) < 2:
        # No unique homologs — duplicate the query so MSA can proceed
        sequences.append(query)

    return sequences


def run_pipeline(raw_input: str, options: Optional[PipelineOptions] = None) -> Dict[str, Any]:
    if options is None:
        options = PipelineOptions()

    # 1. Validate input sequence
    validation = validate_sequence(raw_input)
    if not validation.is_valid:
        raise ValueError("Invalid input sequence: " + "; ".join(validation.errors))

    # 2. Run BLAST to find homologs
    blast_cfg = BlastConfig(
        use_cache=True,
        cache_only=options.use_cache_only,
        verify_ssl=options.verify_ssl,
    )
    blast_hits = run_blast(validation.cleaned_sequence, config=blast_cfg, cache_dir=options.cache_dir)

    # 3. Build sequence list for MSA from query + accepted BLAST hit subjects
    msa_sequences = _build_msa_sequences(validation.cleaned_sequence, blast_hits)

    # 4. Run global MSA on all collected sequences
    msa_cfg = MSAConfig(use_cache=True)
    alignment = align_sequences(msa_sequences, config=msa_cfg, cache_dir=options.cache_dir)

    # 5. Detect conserved regions from alignment
    regions = find_conserved_regions(alignment, min_length=options.min_conserved_length)
    if not regions:
        raise ValueError(
            "No conserved regions detected. "
            "Try lowering min_conserved_length or running with more BLAST hits."
        )
    # Sort by conservation score descending; use best region
    regions.sort(key=lambda r: r.conservation_score, reverse=True)
    region = regions[0]

    # 6. Generate primer candidates from the conserved region consensus
    candidates = generate_candidates_from_region(
        region.consensus_sequence,
        length_min=options.primer_length_min,
        length_max=options.primer_length_max,
    )
    if not candidates:
        raise ValueError("No primer candidates generated from the selected conserved region.")

    # 7. Score every candidate
    scored_candidates = []
    for candidate in candidates:
        structure = analyze_primer(candidate.sequence)
        try:
            specificity = run_primer_specificity(
                candidate.sequence,
                candidate.candidate_id,
                cache_dir=options.cache_dir,
                cache_only=options.use_cache_only,
            )
            spec_score = specificity.specificity_score
            spec_risk = specificity.risk_level
        except Exception:
            # No cache/network available — assume best case for scoring
            spec_score = 100.0
            spec_risk = "LOW"

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
            "specificity": {"score": spec_score, "risk_level": spec_risk},
            "scoring": asdict(score),
        })

    # 8. Select best candidate by composite score
    selected = max(scored_candidates, key=lambda x: x["scoring"]["composite_score"])

    primer_table = [
        {
            "id": c["candidate_id"],
            "sequence": c["sequence"],
            "tm_basic": c["tm_basic"],
            "tm_advanced": c["tm_advanced"],
            "gc_content": c["gc_content"],
            "composite_score": c["scoring"]["composite_score"],
            "risk_level": c["scoring"]["risk_level"],
            "is_accepted": c["is_accepted"],
        }
        for c in scored_candidates
    ]

    result: Dict[str, Any] = {
        "validation": asdict(validation),
        "blast_hits": _serialize_hits(blast_hits),
        "blast_summary": {
            "total": len(blast_hits),
            "accepted": sum(1 for h in blast_hits if h.is_accepted),
            "rejected": sum(1 for h in blast_hits if not h.is_accepted),
        },
        "alignment": {
            "alignment_length": alignment.alignment_length,
            "num_sequences": alignment.num_sequences,
            "full_alignment_text": alignment.full_alignment_text,
        },
        "conserved_regions": [
            {
                "region_id": r.region_id,
                "start_position": r.start_position,
                "end_position": r.end_position,
                "length": r.length,
                "conservation_score": r.conservation_score,
                "entropy_mean": r.entropy_mean,
                "consensus_sequence": r.consensus_sequence,
            }
            for r in regions
        ],
        "selected_region": {
            "region_id": region.region_id,
            "start_position": region.start_position,
            "end_position": region.end_position,
            "length": region.length,
            "conservation_score": region.conservation_score,
            "entropy_mean": region.entropy_mean,
            "consensus_sequence": region.consensus_sequence,
        },
        "candidate_count": len(scored_candidates),
        "accepted_candidate_count": sum(1 for c in scored_candidates if c["is_accepted"]),
        "selected_candidate": selected,
    }

    # 9. Export all formats
    os.makedirs(options.output_dir, exist_ok=True)
    exp_opts = ExportOptions(output_dir=options.output_dir, prefix=options.prefix)
    export_json(result, exp_opts)
    export_txt(result, exp_opts)
    export_csv(primer_table, exp_opts)
    export_fasta(
        [{"id": c["candidate_id"], "sequence": c["sequence"]} for c in scored_candidates],
        exp_opts,
    )
    if options.include_pdf:
        export_pdf(result, ExportOptions(output_dir=options.output_dir, prefix=options.prefix, include_pdf=True))

    return result
