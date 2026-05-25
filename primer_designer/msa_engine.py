import os
import json
import hashlib
import shutil
import subprocess
from dataclasses import dataclass, asdict
from typing import List, Optional

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
        m.update(s.encode('utf-8'))
        m.update(b"\n")
    m.update(json.dumps(asdict(config), sort_keys=True).encode('utf-8'))
    return m.hexdigest()


def _cache_path(cache_dir: str, qhash: str) -> str:
    return os.path.join(cache_dir, f"msa_{qhash}.json")


def align_sequences(sequences: List[str], config: MSAConfig = None, cache_dir: str = "cache") -> AlignmentMatrix:
    """Align sequences using Clustal Omega or ClustalW (global), with optional cache.

    Sequences should be raw sequence strings (no headers). For cache tests, provide
    a precomputed cache file under `cache_dir` named `msa_<hash>.json`.
    """
    if config is None:
        config = MSAConfig()

    _ensure_cache_dir(cache_dir)
    qh = _sequences_hash(sequences, config)
    path = _cache_path(cache_dir, qh)

    if config.use_cache and os.path.exists(path):
        with open(path, 'r', encoding='utf-8') as fh:
            data = json.load(fh)
        return AlignmentMatrix(
            aligned_sequences=data.get('aligned_sequences', []),
            alignment_length=data.get('alignment_length', 0),
            num_sequences=data.get('num_sequences', len(data.get('aligned_sequences', []))),
            full_alignment_text=data.get('full_alignment_text', ''),
        )

    # Try to find aligner
    aligner_exec = None
    if config.aligner == 'clustalo':
        aligner_exec = shutil.which('clustalo') or shutil.which('clustalo.exe')
    if not aligner_exec and config.aligner == 'clustalw' or not aligner_exec:
        # try clustalw2 or clustalw
        aligner_exec = aligner_exec or shutil.which('clustalw2') or shutil.which('clustalw')

    if not aligner_exec:
        raise RuntimeError('No global MSA tool found (clustalo or clustalw). Install Clustal Omega or ClustalW or provide cached alignment.')

    # write sequences to temp fasta
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        fasta_path = os.path.join(td, 'sequences.fasta')
        with open(fasta_path, 'w', encoding='utf-8') as fh:
            for i, s in enumerate(sequences, start=1):
                fh.write(f">seq_{i}\n")
                fh.write(s + "\n")
        out_path = os.path.join(td, 'alignment.aln')
        # build command
        if 'clustalo' in os.path.basename(aligner_exec).lower():
            cmd = [aligner_exec, '--in', fasta_path, '--out', out_path, '--force']
        else:
            # clustalw
            cmd = [aligner_exec, '-INFILE=' + fasta_path, '-OUTFILE=' + out_path, '-OUTPUT=' + config.output_format.upper(), '-ALIGN']
        subprocess.run(cmd, check=True)
        with open(out_path, 'r', encoding='utf-8') as fh:
            aln_text = fh.read()
        # parse simple aligned sequences from clustal output: collect lines without headers
        aligned = []
        for line in aln_text.splitlines():
            if line and not line.startswith(' ') and '\t' not in line and not line.startswith('CLUSTAL'):
                parts = line.split()
                if len(parts) >= 2 and parts[0].startswith('seq_'):
                    aligned.append(parts[1])
        alignment_length = len(aligned[0]) if aligned else 0

    # save to cache
    tosave = {
        'aligned_sequences': aligned,
        'alignment_length': alignment_length,
        'num_sequences': len(aligned),
        'full_alignment_text': aln_text,
    }
    with open(path, 'w', encoding='utf-8') as fh:
        json.dump(tosave, fh, indent=2)

    return AlignmentMatrix(aligned, alignment_length, len(aligned), aln_text)
