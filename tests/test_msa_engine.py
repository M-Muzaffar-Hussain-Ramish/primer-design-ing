import json
from primer_designer.msa_engine import align_sequences, MSAConfig, _sequences_hash, _cache_path


def test_load_alignment_from_cache(tmp_path):
    seqs = ["ATGCGTACGTTAGC", "ATGCGTACGTTAGC"]
    cfg = MSAConfig(use_cache=True)
    qh = _sequences_hash(seqs, cfg)
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    path = cache_dir / _cache_path(str(cache_dir), qh).split('/')[-1]

    sample = {
        'aligned_sequences': ['ATGCGTACGTTAGC', 'ATGCGTACGTTAGC'],
        'alignment_length': 13,
        'num_sequences': 2,
        'full_alignment_text': 'CLUSTAL W (1.83) multiple sequence alignment\n\nseq_1    ATGCGTACGTTAGC\nseq_2    ATGCGTACGTTAGC\n'
    }
    with open(path, 'w', encoding='utf-8') as fh:
        json.dump(sample, fh)

    mat = align_sequences(seqs, config=cfg, cache_dir=str(cache_dir))
    assert mat.num_sequences == 2
    assert mat.alignment_length == 13
    assert 'CLUSTAL' in mat.full_alignment_text


def test_biopython_fallback_when_no_external_aligner(tmp_path):
    """When no external aligner (clustalo/clustalw) is on PATH, BioPython star-alignment is used."""
    seqs = ["ATGCGTACGT", "ATGCGTACGT"]
    cfg = MSAConfig(use_cache=False, aligner='clustalo')
    import os
    old_path = os.environ.get('PATH')
    os.environ['PATH'] = ''
    try:
        mat = align_sequences(seqs, config=cfg, cache_dir=str(tmp_path))
        assert mat.num_sequences == 2
        assert mat.alignment_length > 0
        assert "BioPython" in mat.full_alignment_text
    finally:
        os.environ['PATH'] = old_path
