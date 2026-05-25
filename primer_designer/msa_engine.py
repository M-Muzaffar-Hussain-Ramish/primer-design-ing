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
    """Parse Clustal/Clustal Omega output, handling multi-block alignments.

    Multi-block format repeats sequence names every ~60 columns; this
    accumulates all blocks per sequence name and returns sequences in
    original order.
    """
    seq_dict: Dict[str, str] = {}
    seq_order: List[str] = []

    for line in aln_text.splitlines():
        line = line.rstrip()
        # Skip blank lines, CLUSTAL header, and conservation lines (start with space)
        if not line or line.startswith("CLUSTAL") or line.startswith(" "):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        name = parts[0]
        # The sequence part is everything after the name; strip trailing position numbers
        seq_part = parts[1]
        # Remove any trailing digit-only token (position counter)
        if len(parts) >= 3 and parts[-1].isdigit():
            seq_part = "".join(parts[1:-1])
        else:
            seq_part = "".join(parts[1:])
        # Only keep valid alignment characters
        seq_part = "".join(c for c in seq_part if c.isalpha() or c == "-")
        if not seq_part:
            continue
        if name not in seq_dict:
            seq_dict[name] = ""
            seq_order.append(name)
        seq_dict[name] += seq_part

    return [seq_dict[n] for n in seq_order]


def align_sequences(sequences: List[str], config: MSAConfig = None, cache_dir: str = "cache") -> AlignmentMatrix:
    """Align sequences using Clustal Omega or ClustalW (global MSA), with optional cache.

    Pass raw sequence strings without headers. Pre-populate the cache directory
    with a ``msa_<hash>.json`` file to skip the external tool call (useful in tests/CI).

    Environment:
        CLUSTALO / CLUSTALW env vars can override binary discovery.
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

    # Discover aligner binary
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

    if not aligner_exec:
        raise RuntimeError(
            "No global MSA aligner found (clustalo or clustalw). "
            "Install Clustal Omega: https://www.ebi.ac.uk/Tools/msa/clustalo/ "
            "or ClustalW: http://www.clustal.org/clustal2/. "
            "Alternatively supply a pre-computed cache file."
        )

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
        except subprocess.CalledProcessError as exc:
            stderr = exc.stderr.decode(errors="replace") if exc.stderr else ""
            raise RuntimeError(f"MSA aligner failed: {stderr}") from exc

        with open(out_path, "r", encoding="utf-8") as fh:
            aln_text = fh.read()

    aligned = _parse_clustal(aln_text)
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
