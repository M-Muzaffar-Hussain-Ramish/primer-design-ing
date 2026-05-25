import os
import json
import hashlib
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional


@dataclass
class MSAConfig:
    aligner: str = "clustalo"
    output_format: str = "clustal"
    use_cache: bool = True


@dataclass
class AlignmentMatrix:
    aligned_sequences: List[str]
    alignment_length: int
    num_sequences: int
    full_alignment_text: str


def _ensure_cache_dir(cache_dir: str):
    os.makedirs(cache_dir, exist_ok=True)


def _sequences_hash(sequences: List[str], config: MSAConfig) -> str:
    m = hashlib.sha256()
    for s in sequences:
        m.update(s.encode("utf-8"))
        m.update(b"\n")
    m.update(json.dumps(asdict(config), sort_keys=True).encode("utf-8"))
    return m.hexdigest()


def _cache_path(cache_dir: str, qhash: str) -> str:
    return os.path.join(cache_dir, f"msa_{qhash}.json")


def _parse_clustal(aln_text: str) -> List[str]:
    """Parse Clustal/Clustal Omega output handling multi-block alignments."""
    seq_dict: Dict[str, str] = {}
    seq_order: List[str] = []

    for line in aln_text.splitlines():
        line = line.rstrip()
        if not line or line.startswith("CLUSTAL") or line.startswith(" "):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        name = parts[0]
        seq_part = parts[1]
        if len(parts) >= 3 and parts[-1].isdigit():
            seq_part = "".join(parts[1:-1])
        else:
            seq_part = "".join(parts[1:])
        seq_part = "".join(c for c in seq_part if c.isalpha() or c == "-")
        if not seq_part:
            continue
        if name not in seq_dict:
            seq_dict[name] = ""
            seq_order.append(name)
        seq_dict[name] += seq_part

    return [seq_dict[n] for n in seq_order]


def _biopython_star_msa(sequences: List[str]) -> List[str]:
    """Star MSA via BioPython pairwise aligner — no external tool required.

    Aligns every sequence against the first (query) sequence as reference.
    Gap characters are inserted to produce equal-length aligned strings.
    Falls back gracefully if BioPython is not available.
    """
    try:
        from Bio import pairwise2  # type: ignore
        _USE_PW2 = True
    except ImportError:
        _USE_PW2 = False

    if not _USE_PW2:
        try:
            from Bio.Align import PairwiseAligner  # type: ignore
            _USE_PALIGN = True
        except ImportError:
            _USE_PALIGN = False
    else:
        _USE_PALIGN = False

    if not _USE_PW2 and not _USE_PALIGN:
        raise RuntimeError(
            "BioPython is required for MSA fallback. Install with: pip install biopython"
        )

    ref = sequences[0]
    aligned: List[str] = []

    for seq in sequences:
        if seq == ref:
            aligned.append(ref)
            continue

        if _USE_PW2:
            from Bio import pairwise2  # type: ignore
            alns = pairwise2.align.globalms(ref, seq, 2, -1, -3, -0.5)
            if alns:
                aligned.append(str(alns[0].seqB))
            else:
                aligned.append(seq)
        else:
            from Bio.Align import PairwiseAligner  # type: ignore
            aligner = PairwiseAligner()
            aligner.mode = "global"
            aligner.match_score = 2
            aligner.mismatch_score = -1
            aligner.open_gap_score = -3
            aligner.extend_gap_score = -0.5
            alns = list(aligner.align(ref, seq))
            if not alns:
                aligned.append(seq)
                continue
            best = alns[0]
            # Build gap-annotated string from alignment coordinates
            aligned.append(_extract_query_aligned(ref, seq, best.aligned))

    return aligned


def _extract_query_aligned(ref: str, seq: str, aligned_coords) -> str:
    """Build the gap-annotated query string from BioPython PairwiseAligner coordinates."""
    ref_blocks = list(aligned_coords[0])   # (start, end) in ref
    seq_blocks = list(aligned_coords[1])   # (start, end) in seq

    result: List[str] = []
    r_pos = 0
    s_pos = 0

    for (r_start, r_end), (s_start, s_end) in zip(ref_blocks, seq_blocks):
        # Gap in seq (deletion relative to ref)
        if r_start > r_pos:
            result.extend(["-"] * (r_start - r_pos))
        # Gap in ref (insertion relative to ref) — include extra seq bases
        if s_start > s_pos:
            result.extend(list(seq[s_pos:s_start]))
        # Aligned block
        result.extend(list(seq[s_start:s_end]))
        r_pos = r_end
        s_pos = s_end

    # Trailing
    if r_pos < len(ref):
        result.extend(["-"] * (len(ref) - r_pos))
    if s_pos < len(seq):
        result.extend(list(seq[s_pos:]))

    return "".join(result)


def align_sequences(sequences: List[str], config: MSAConfig = None, cache_dir: str = "cache") -> AlignmentMatrix:
    """Align sequences using Clustal Omega/ClustalW (preferred) or BioPython fallback.

    Pass raw DNA strings without FASTA headers.
    Pre-populate ``cache_dir/msa_<hash>.json`` to skip alignment entirely (useful in tests/CI).

    Discovery order:
      1. CLUSTALO or CLUSTALW environment variable (path to binary)
      2. clustalo / clustalo.exe / clustalw2 / clustalw on PATH
      3. BioPython pairwise star-alignment fallback (no external tool needed)

    Set NCBI_VERIFY_SSL=0 in your .env for corporate proxy environments.
    """
    if config is None:
        config = MSAConfig()

    _ensure_cache_dir(cache_dir)
    qh = _sequences_hash(sequences, config)
    path = _cache_path(cache_dir, qh)

    if config.use_cache and os.path.exists(path):
        with open(path, "r", encoding="utf-8") as fh:
            data = json.load(fh)
        return AlignmentMatrix(
            aligned_sequences=data.get("aligned_sequences", []),
            alignment_length=data.get("alignment_length", 0),
            num_sequences=data.get("num_sequences", len(data.get("aligned_sequences", []))),
            full_alignment_text=data.get("full_alignment_text", ""),
        )

    # ── Try external aligner ──────────────────────────────────────────────────
    aligner_exec: Optional[str] = None
    env_override = os.environ.get("CLUSTALO") or os.environ.get("CLUSTALW")
    if env_override and os.path.isfile(env_override):
        aligner_exec = env_override
    else:
        aligner_exec = (
            shutil.which("clustalo")
            or shutil.which("clustalo.exe")
            or shutil.which("clustalw2")
            or shutil.which("clustalw")
            or shutil.which("clustalw.exe")
        )

    aligned: List[str] = []
    aln_text: str = ""

    if aligner_exec:
        with tempfile.TemporaryDirectory() as td:
            fasta_path = os.path.join(td, "sequences.fasta")
            with open(fasta_path, "w", encoding="utf-8") as fh:
                for i, s in enumerate(sequences, start=1):
                    fh.write(f">seq_{i}\n{s}\n")
            out_path = os.path.join(td, "alignment.aln")

            basename = os.path.basename(aligner_exec).lower()
            if "clustalo" in basename:
                cmd = [aligner_exec, "--in", fasta_path, "--out", out_path, "--force", "--outfmt=clustal"]
            else:
                cmd = [aligner_exec, f"-INFILE={fasta_path}", f"-OUTFILE={out_path}", "-OUTPUT=CLUSTAL", "-ALIGN"]

            try:
                subprocess.run(cmd, check=True, capture_output=True)
                with open(out_path, "r", encoding="utf-8") as fh:
                    aln_text = fh.read()
                aligned = _parse_clustal(aln_text)
            except (subprocess.CalledProcessError, OSError):
                aligner_exec = None  # fall through to BioPython

    if not aligned:
        # ── BioPython pairwise fallback (works on any OS without external tools) ─
        try:
            aligned = _biopython_star_msa(sequences)
            aln_text = (
                "BioPython star alignment (pairwise2 fallback — install clustalo for true MSA)\n\n"
                + "\n".join(f"seq_{i+1}    {s}" for i, s in enumerate(aligned))
            )
        except Exception as exc:
            raise RuntimeError(
                "MSA failed: no external aligner (clustalo/clustalw) found and BioPython fallback "
                "also failed. Install Clustal Omega: https://www.ebi.ac.uk/Tools/msa/clustalo/ "
                f"Error: {exc}"
            ) from exc

    alignment_length = max((len(s) for s in aligned), default=0)

    tosave = {
        "aligned_sequences": aligned,
        "alignment_length": alignment_length,
        "num_sequences": len(aligned),
        "full_alignment_text": aln_text,
    }
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(tosave, fh, indent=2)

    return AlignmentMatrix(aligned, alignment_length, len(aligned), aln_text)
