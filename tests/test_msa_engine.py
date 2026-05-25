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


def test_missing_aligner_raises(tmp_path):
    seqs = ["ATGC", "ATGC"]
    cfg = MSAConfig(use_cache=False, aligner='clustalo')
    # ensure no clustalo/clustalw in PATH by using a temp PATH
    import os
    old_path = os.environ.get('PATH')
    os.environ['PATH'] = ''
    try:
        raised = False
        try:
            align_sequences(seqs, config=cfg, cache_dir=str(tmp_path))
        except RuntimeError as e:
            raised = True
            assert 'No global MSA' in str(e)
        assert raised
    finally:
        os.environ['PATH'] = old_path
