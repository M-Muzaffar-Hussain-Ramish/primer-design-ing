import json
import os
from dataclasses import asdict, dataclass
from typing import Any, Dict, List, Optional, Tuple

from primer_designer.blast_engine import BlastConfig, BlastHit, run_blast
from primer_designer.conserved_region_detector import find_conserved_regions
from primer_designer.export_manager import ExportOptions, export_csv, export_fasta, export_json, export_pdf, export_txt
from primer_designer.msa_engine import MSAConfig, align_sequences
from primer_designer.primer_generator import (
    generate_candidates_from_region,
    generate_reverse_candidates_from_region,
)
from primer_designer.scoring_engine import score_primer
from primer_designer.secondary_structure import analyze_primer
from primer_designer.sequence_validator import validate_sequence
from primer_designer.specificity_checker import run_primer_specificity
from primer_designer.thermodynamics import thermo_profile


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
    blast_max_hits: int = 20           # max homolog sequences to retrieve
    blast_identity_min: float = 95.0   # minimum identity % for BLAST hits used in MSA
    max_tm_diff: float = 2.0           # max allowed Tm difference between F and R primers


def _serialize_hits(hits: List[BlastHit]) -> List[Dict[str, Any]]:
    return [asdict(hit) for hit in hits]


def _build_msa_sequences(query: str, blast_hits: List[BlastHit], max_seqs: int = 20) -> List[str]:
    """Build sequence list for MSA: query + unique accepted BLAST hit subjects."""
    sequences: List[str] = [query]
    seen = {query}

    for hit in blast_hits:
        if not hit.is_accepted:
            continue
        if len(sequences) >= max_seqs:
            break
        subj_raw = hit.subject_sequence or hit.aligned_segment_preview
        subj_clean = "".join(c for c in subj_raw if c.isalpha())
        if len(subj_clean) < 18 or subj_clean in seen:
            continue
        sequences.append(subj_clean)
        seen.add(subj_clean)

    if len(sequences) < 2:
        sequences.append(query)  # duplicate query so MSA has ≥2 sequences

    return sequences


def _score_candidate(candidate, conservation_score: float, specificity_score: float, risk: str):
    """Score a primer candidate and attach structure analysis."""
    structure = analyze_primer(candidate.sequence)
    score = score_primer(
        primer_id=candidate.candidate_id,
        conservation_score=conservation_score,
        tm_basic=candidate.tm_basic,
        tm_advanced=candidate.tm_advanced,
        specificity_score=specificity_score,
        secondary_structure_risk=structure.overall_risk,
        gc_content=candidate.gc_content,
        length=candidate.length,
    )
    return structure, score


def _select_best_pair(
    fwd_candidates: List, rev_candidates: List, conservation_score: float, max_tm_diff: float
) -> Tuple[Optional[Dict], Optional[Dict]]:
    """Select the best accepted F + R primer pair with compatible Tm values.

    Strategy:
      1. Score every accepted forward and reverse candidate independently.
      2. Pick the top-20 forward and top-20 reverse by composite score.
      3. Among those, find the pair with the smallest |Tm_F - Tm_R| that is
         within max_tm_diff, with no cross-dimer fail.
      4. Fall back to best individual candidates if no compatible pair found.
    """
    def enrich(candidates, spec_score=100.0, spec_risk="LOW"):
        enriched = []
        for c in candidates:
            if not c.is_accepted:
                continue
            st, sc = _score_candidate(c, conservation_score, spec_score, spec_risk)
            enriched.append({
                **asdict(c),
                "structure": asdict(st),
                "specificity": {"score": spec_score, "risk_level": spec_risk},
                "scoring": asdict(sc),
            })
        enriched.sort(key=lambda x: x["scoring"]["composite_score"], reverse=True)
        return enriched

    top_fwd = enrich(fwd_candidates)[:20]
    top_rev = enrich(rev_candidates)[:20]

    if not top_fwd or not top_rev:
        best_f = top_fwd[0] if top_fwd else None
        best_r = top_rev[0] if top_rev else None
        return best_f, best_r

    # Find best compatible pair
    best_pair = None
    best_pair_score = -1.0

    for f in top_fwd:
        for r in top_rev:
            tm_diff = abs(f["tm_advanced"] - r["tm_advanced"])
            if tm_diff > max_tm_diff:
                continue
            # Check cross-dimer
            cross = analyze_primer(f["sequence"], partner_seq=r["sequence"])
            if cross.pass_fail is False:
                continue
            pair_score = f["scoring"]["composite_score"] + r["scoring"]["composite_score"]
            if pair_score > best_pair_score:
                best_pair_score = pair_score
                best_pair = (f, r)

    if best_pair:
        return best_pair

    # No compatible pair — return best of each independently
    return top_fwd[0], top_rev[0]


def run_pipeline(raw_input: str, options: Optional[PipelineOptions] = None) -> Dict[str, Any]:
    """Full primer design pipeline.

    Flow:
      1. Validate input sequence
      2. BLAST → select homologs ≥ blast_identity_min %
      3. MSA (clustalo preferred, BioPython pairwise fallback)
      4. Detect conserved regions
      5. Generate forward + reverse primer candidates
      6. Score and select best compatible F+R pair (Tm diff ≤ max_tm_diff)
      7. Export results (JSON / TXT / CSV / FASTA / PDF)
    """
    if options is None:
        options = PipelineOptions()

    # ── 1. Validate ────────────────────────────────────────────────────────────
    validation = validate_sequence(raw_input)
    if not validation.is_valid:
        raise ValueError("Invalid input sequence: " + "; ".join(validation.errors))

    # ── 2. BLAST ───────────────────────────────────────────────────────────────
    blast_cfg = BlastConfig(
        identity_threshold=options.blast_identity_min,
        coverage_threshold=70.0,
        max_hits=options.blast_max_hits,
        use_cache=True,
        cache_only=options.use_cache_only,
        verify_ssl=options.verify_ssl,
    )
    blast_hits = run_blast(validation.cleaned_sequence, config=blast_cfg, cache_dir=options.cache_dir)

    accepted_hits = [h for h in blast_hits if h.is_accepted]

    # ── 3. MSA ─────────────────────────────────────────────────────────────────
    msa_sequences = _build_msa_sequences(validation.cleaned_sequence, blast_hits)
    msa_cfg = MSAConfig(use_cache=True)
    alignment = align_sequences(msa_sequences, config=msa_cfg, cache_dir=options.cache_dir)

    # ── 4. Conserved regions ───────────────────────────────────────────────────
    regions = find_conserved_regions(alignment, min_length=options.min_conserved_length)
    if not regions:
        raise ValueError(
            "No conserved regions detected in the MSA. "
            "Try a longer input sequence, lower min_conserved_length, or check BLAST hits."
        )
    regions.sort(key=lambda r: r.conservation_score, reverse=True)
    region = regions[0]  # best conserved region

    # ── 5. Generate F + R candidates ──────────────────────────────────────────
    fwd_candidates = generate_candidates_from_region(
        region.consensus_sequence,
        length_min=options.primer_length_min,
        length_max=options.primer_length_max,
    )
    rev_candidates = generate_reverse_candidates_from_region(
        region.consensus_sequence,
        length_min=options.primer_length_min,
        length_max=options.primer_length_max,
    )

    if not fwd_candidates and not rev_candidates:
        raise ValueError("No primer candidates generated from the selected conserved region.")

    # ── 6. Score + select best F+R pair ───────────────────────────────────────
    selected_forward, selected_reverse = _select_best_pair(
        fwd_candidates, rev_candidates, region.conservation_score, options.max_tm_diff
    )

    if not selected_forward or not selected_reverse:
        raise ValueError(
            "Could not find a valid forward+reverse primer pair. "
            "Try relaxing Tm compatibility or primer length constraints."
        )

    # Amplicon size (distance from start of F to end of R binding site)
    fwd_len = selected_forward["length"]
    rev_len = selected_reverse["length"]
    region_len = region.length
    amplicon_size = region_len  # F binds 5' end, R binds 3' end of conserved region

    tm_f = selected_forward["tm_advanced"]
    tm_r = selected_reverse["tm_advanced"]

    # All scored candidates for the primer table
    def _all_scored(candidates, spec_score=100.0):
        out = []
        for c in candidates:
            st = analyze_primer(c.sequence)
            sc = score_primer(
                primer_id=c.candidate_id,
                conservation_score=region.conservation_score,
                tm_basic=c.tm_basic, tm_advanced=c.tm_advanced,
                specificity_score=spec_score,
                secondary_structure_risk=st.overall_risk,
                gc_content=c.gc_content, length=c.length,
            )
            out.append({**asdict(c), "structure": asdict(st), "scoring": asdict(sc)})
        return out

    all_fwd_scored = _all_scored(fwd_candidates)
    all_rev_scored = _all_scored(rev_candidates)

    primer_table = [
        {
            "id": c["candidate_id"],
            "direction": c["direction"],
            "sequence": c["sequence"],
            "length": c["length"],
            "tm_basic": c["tm_basic"],
            "tm_advanced": c["tm_advanced"],
            "gc_content": c["gc_content"],
            "composite_score": c["scoring"]["composite_score"],
            "risk_level": c["scoring"]["risk_level"],
            "is_accepted": c["is_accepted"],
        }
        for c in all_fwd_scored + all_rev_scored
    ]

    # ── 7. Build result dict ────────────────────────────────────────────────────
    result: Dict[str, Any] = {
        "validation": asdict(validation),
        "blast_summary": {
            "total": len(blast_hits),
            "accepted": len(accepted_hits),
            "rejected": len(blast_hits) - len(accepted_hits),
            "identity_threshold": options.blast_identity_min,
            "sequences_used_for_msa": len(msa_sequences),
        },
        "blast_hits": _serialize_hits(blast_hits),
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
        "primer_pair": {
            "forward": selected_forward,
            "reverse": selected_reverse,
            "tm_forward": round(tm_f, 2),
            "tm_reverse": round(tm_r, 2),
            "tm_difference": round(abs(tm_f - tm_r), 2),
            "tm_compatible": abs(tm_f - tm_r) <= options.max_tm_diff,
            "amplicon_size_bp": amplicon_size,
        },
        "forward_candidate_count": len(fwd_candidates),
        "reverse_candidate_count": len(rev_candidates),
        "accepted_forward_count": sum(1 for c in fwd_candidates if c.is_accepted),
        "accepted_reverse_count": sum(1 for c in rev_candidates if c.is_accepted),
        # Legacy key kept for backward compatibility
        "selected_candidate": selected_forward,
        "candidate_count": len(fwd_candidates) + len(rev_candidates),
        "accepted_candidate_count": sum(1 for c in fwd_candidates + rev_candidates if c.is_accepted),
    }

    # ── 8. Export ───────────────────────────────────────────────────────────────
    os.makedirs(options.output_dir, exist_ok=True)
    exp_opts = ExportOptions(output_dir=options.output_dir, prefix=options.prefix)
    export_json(result, exp_opts)
    export_txt(result, exp_opts)
    export_csv(primer_table, exp_opts)
    export_fasta(
        [
            {"id": selected_forward["candidate_id"], "sequence": selected_forward["sequence"]},
            {"id": selected_reverse["candidate_id"], "sequence": selected_reverse["sequence"]},
        ],
        exp_opts,
    )
    if options.include_pdf:
        export_pdf(result, ExportOptions(output_dir=options.output_dir, prefix=options.prefix, include_pdf=True))

    return result
