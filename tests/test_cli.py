"""Tests for the command-line interface."""
import subprocess
import sys
from pathlib import Path


PYTHON = sys.executable
SEQ = "ATGCGTACGTTAGCCTAGCTATGCGTACGTTAGCCTAGCT"


def _run(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [PYTHON, "-m", "primer_designer.cli", *args],
        capture_output=True,
        text=True,
    )


def test_cli_raw_sequence():
    result = _run(SEQ)
    assert result.returncode == 0
    assert "Forward" in result.stdout or "Reverse" in result.stdout


def test_cli_fasta_file(tmp_path: Path):
    fasta = tmp_path / "test.fa"
    fasta.write_text(f">test_seq\n{SEQ}\n")
    result = _run(str(fasta))
    assert result.returncode == 0
    assert "Forward" in result.stdout or "Reverse" in result.stdout


def test_cli_custom_length():
    result = _run(SEQ, "--length", "18")
    assert result.returncode == 0


def test_cli_invalid_sequence():
    result = _run("ATGCXYZ!!!")
    # Should either fail gracefully (non-zero exit or stderr) or print validation error
    output = result.stdout + result.stderr
    assert "invalid" in output.lower() or "error" in output.lower() or result.returncode != 0


def test_cli_too_short():
    result = _run("ATGC")
    output = result.stdout + result.stderr
    assert "short" in output.lower() or "error" in output.lower() or result.returncode != 0
