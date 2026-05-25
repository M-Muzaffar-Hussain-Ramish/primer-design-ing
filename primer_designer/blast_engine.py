import json
import os
import hashlib
from dataclasses import dataclass, asdict
from typing import List, Optional

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


def _ensure_cache_dir(cache_dir: str):
    os.makedirs(cache_dir, exist_ok=True)


def _query_hash(query: str, config: BlastConfig) -> str:
    m = hashlib.sha256()
    m.update(query.encode("utf-8"))
    # include config values relevant to output
    m.update(json.dumps(asdict(config), sort_keys=True).encode("utf-8"))
    return m.hexdigest()


def _cache_path(cache_dir: str, qhash: str) -> str:
    return os.path.join(cache_dir, f"{qhash}.json")


def run_blast(query: str, config: BlastConfig = None, cache_dir: str = "cache") -> List[BlastHit]:
    """Run BLAST (or load from cache). When `config.use_cache` is True and a cached
    JSON exists for the query+config, it will be loaded and returned.

    NOTE: This function will attempt an online BLASTn if no cache is available and
    Biopython is installed, but network access may fail in tests/CI. Prefer using
    `use_cache=True` with prepopulated cache during testing.
    """
    if config is None:
        config = BlastConfig()

    _ensure_cache_dir(cache_dir)
    qh = _query_hash(query, config)
    path = _cache_path(cache_dir, qh)

    if config.use_cache and os.path.exists(path):
        with open(path, "r", encoding="utf-8") as fh:
            data = json.load(fh)
        hits = [BlastHit(**h) for h in data.get("hits", [])]
        return hits

    if config.use_cache and config.cache_only:
        raise RuntimeError(f"Cache-only BLAST requested but no cache entry found for query hash {qh}")

    # No cache: try to perform online BLAST using Biopython
    try:
        from Bio.Blast import NCBIWWW
        from Bio.Blast import NCBIXML
    except Exception:
        raise RuntimeError("Biopython not available or BLAST runtime not configured; provide cached results")

    # run qblast
    res = NCBIWWW.qblast(config.program, config.database, query)
    blast_records = NCBIXML.parse(res)

    hits: List[BlastHit] = []
    for record in blast_records:
        for aln in record.alignments:
            for hsp in aln.hsps:
                # basic parsing; break if too many
                hit = BlastHit(
                    accession=getattr(aln, 'accession', getattr(aln, 'hit_id', '')),
                    organism=aln.hit_def if hasattr(aln, 'hit_def') else "",
                    alignment_length=hsp.align_length,
                    query_coverage=100.0 * (hsp.align_length / len(query)) if len(query) else 0.0,
                    identity_percent=100.0 * (hsp.identities / hsp.align_length) if hsp.align_length else 0.0,
                    mismatches=hsp.align_length - hsp.identities,
                    gaps=getattr(hsp, 'gaps', 0) or 0,
                    e_value=float(hsp.expect),
                    bit_score=float(hsp.bits),
                    aligned_segment_preview=(hsp.query[:80] if hasattr(hsp, 'query') else ""),
                    query_start=getattr(hsp, 'query_start', 0),
                    query_end=getattr(hsp, 'query_end', 0),
                    subject_start=getattr(hsp, 'sbjct_start', 0),
                    subject_end=getattr(hsp, 'sbjct_end', 0),
                    strand="plus/plus",
                    is_accepted=True,
                    rejection_reason=None,
                    raw_score=float(getattr(hsp, 'score', 0.0)),
                )
                hits.append(hit)
                if len(hits) >= config.max_hits:
                    break
            if len(hits) >= config.max_hits:
                break
        if len(hits) >= config.max_hits:
            break

    # save to cache
    tosave = {"query_hash": qh, "hits": [asdict(h) for h in hits]}
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(tosave, fh, indent=2)

    return hits
