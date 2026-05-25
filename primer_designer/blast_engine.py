import json
import os
import hashlib
import time
from dataclasses import dataclass, asdict, field
from typing import List, Optional

try:
    from loguru import logger
except ImportError:
    import logging as _logging
    logger = _logging.getLogger(__name__)


@dataclass
class BlastConfig:
    database: str = "nt"
    program: str = "blastn"
    identity_threshold: float = 95.0
    coverage_threshold: float = 80.0
    max_hits: int = 100
    e_value_cutoff: float = 1e-10
    entrez_query: str = ""
    use_cache: bool = True
    cache_only: bool = False
    verify_ssl: bool = True


@dataclass
class BlastHit:
    accession: str
    organism: str
    alignment_length: int
    query_coverage: float
    identity_percent: float
    mismatches: int
    gaps: int
    e_value: float
    bit_score: float
    aligned_segment_preview: str
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    strand: str
    is_accepted: bool
    rejection_reason: Optional[str]
    raw_score: float
    subject_sequence: str = ""  # full aligned subject sequence from BLAST


def _ensure_cache_dir(cache_dir: str):
    os.makedirs(cache_dir, exist_ok=True)


def _query_hash(query: str, config: BlastConfig) -> str:
    m = hashlib.sha256()
    m.update(query.encode("utf-8"))
    m.update(json.dumps(asdict(config), sort_keys=True).encode("utf-8"))
    return m.hexdigest()


def _cache_path(cache_dir: str, qhash: str) -> str:
    return os.path.join(cache_dir, f"{qhash}.json")


def _apply_filters(hits: List[BlastHit], config: BlastConfig) -> List[BlastHit]:
    """Tag each hit as accepted or rejected based on config thresholds."""
    for hit in hits:
        reasons = []
        if hit.identity_percent < config.identity_threshold:
            reasons.append(f"identity {hit.identity_percent:.1f}% < {config.identity_threshold}%")
        if hit.query_coverage < config.coverage_threshold:
            reasons.append(f"coverage {hit.query_coverage:.1f}% < {config.coverage_threshold}%")
        if hit.e_value > config.e_value_cutoff:
            reasons.append(f"e_value {hit.e_value:.2e} > {config.e_value_cutoff:.2e}")
        if reasons:
            hit.is_accepted = False
            hit.rejection_reason = "; ".join(reasons)
        else:
            hit.is_accepted = True
            hit.rejection_reason = None
    return hits


def _run_qblast_with_retry(program: str, database: str, query: str, entrez_query: str = "", max_retries: int = 3):
    """Call NCBIWWW.qblast with exponential backoff on transient errors."""
    from Bio.Blast import NCBIWWW

    api_key = os.environ.get("NCBI_API_KEY", "")
    kwargs = {}
    if api_key:
        kwargs["entrez_query"] = entrez_query
        # Biopython passes ncbi_gi, format_type etc; pass api_key via hitlist_size workaround not needed
        # qblast doesn't expose api_key directly but we set it via env for some versions
        os.environ["NCBI_API_KEY"] = api_key

    delays = [5, 30, 120]
    last_exc = None
    for attempt in range(max_retries):
        try:
            result = NCBIWWW.qblast(program, database, query, entrez_query=entrez_query)
            return result
        except Exception as exc:
            last_exc = exc
            logger.warning(f"BLAST attempt {attempt + 1}/{max_retries} failed: {exc}")
            if attempt < max_retries - 1:
                time.sleep(delays[attempt])
    raise RuntimeError(f"BLAST failed after {max_retries} attempts: {last_exc}") from last_exc


def run_blast(query: str, config: BlastConfig = None, cache_dir: str = "cache") -> List[BlastHit]:
    """Run BLAST (or load from cache). Hits are tagged is_accepted/rejection_reason per config thresholds.

    Set NCBI_API_KEY env var for higher rate limits (10 req/s vs 3 req/s).
    Set config.cache_only=True to raise RuntimeError instead of making network calls.
    """
    if config is None:
        config = BlastConfig()

    _ensure_cache_dir(cache_dir)
    qh = _query_hash(query, config)
    path = _cache_path(cache_dir, qh)

    if config.use_cache and os.path.exists(path):
        with open(path, "r", encoding="utf-8") as fh:
            data = json.load(fh)
        hits = []
        for h in data.get("hits", []):
            # Backcompat: old caches may not have subject_sequence
            h.setdefault("subject_sequence", "")
            hits.append(BlastHit(**h))
        # Re-apply filters so cached hits reflect current config thresholds
        hits = _apply_filters(hits, config)
        return hits

    if config.cache_only:
        raise RuntimeError(f"Cache-only BLAST requested but no cache entry found for query hash {qh}")

    ssl_module = None
    original_context = None
    if not config.verify_ssl:
        try:
            import ssl as ssl_module
            original_context = ssl_module._create_default_https_context
            ssl_module._create_default_https_context = ssl_module._create_unverified_context
            os.environ["PYTHONHTTPSVERIFY"] = "0"
        except Exception:
            ssl_module = None

    try:
        from Bio.Blast import NCBIXML
    except Exception as exc:
        raise RuntimeError(
            "Biopython is required for live BLAST. Install with: pip install biopython. "
            "Alternatively supply a pre-computed cache file."
        ) from exc

    try:
        res = _run_qblast_with_retry(config.program, config.database, query, config.entrez_query)
        blast_records = NCBIXML.parse(res)
    finally:
        if ssl_module and original_context is not None:
            ssl_module._create_default_https_context = original_context
            os.environ.pop("PYTHONHTTPSVERIFY", None)

    hits: List[BlastHit] = []
    for record in blast_records:
        for aln in record.alignments:
            for hsp in aln.hsps:
                query_len = len(query) or 1
                coverage = 100.0 * (hsp.align_length / query_len)
                identity = 100.0 * (hsp.identities / hsp.align_length) if hsp.align_length else 0.0
                hit = BlastHit(
                    accession=getattr(aln, "accession", getattr(aln, "hit_id", "")),
                    organism=aln.hit_def if hasattr(aln, "hit_def") else "",
                    alignment_length=hsp.align_length,
                    query_coverage=coverage,
                    identity_percent=identity,
                    mismatches=hsp.align_length - hsp.identities,
                    gaps=getattr(hsp, "gaps", 0) or 0,
                    e_value=float(hsp.expect),
                    bit_score=float(hsp.bits),
                    aligned_segment_preview=(hsp.query[:80] if hasattr(hsp, "query") else ""),
                    query_start=getattr(hsp, "query_start", 0),
                    query_end=getattr(hsp, "query_end", 0),
                    subject_start=getattr(hsp, "sbjct_start", 0),
                    subject_end=getattr(hsp, "sbjct_end", 0),
                    strand="plus/plus",
                    is_accepted=True,
                    rejection_reason=None,
                    raw_score=float(getattr(hsp, "score", 0.0)),
                    subject_sequence=getattr(hsp, "sbjct", "").replace("-", "") if hasattr(hsp, "sbjct") else "",
                )
                hits.append(hit)
                if len(hits) >= config.max_hits:
                    break
            if len(hits) >= config.max_hits:
                break
        if len(hits) >= config.max_hits:
            break

    hits = _apply_filters(hits, config)

    # Save to cache
    tosave = {"query_hash": qh, "hits": [asdict(h) for h in hits]}
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(tosave, fh, indent=2)

    return hits
