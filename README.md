Primer Designer - Minimal Implementation

This repository contains a minimal, well-tested subset of the Research-Grade PCR Primer Design System described in DESIGN.md.

Quickstart:

1. Create a virtualenv and install test dependency:

```bash
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
```

2. Run tests:

```bash
pytest -q
```

3. Run CLI on a raw sequence:

```bash
python -m primer_designer.cli "ATGCGTACGTTAGCCTAGCT..."
```
