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

4. Run the Streamlit UI:

```bash
streamlit run streamlit_app.py
```

5. If you prefer the app folder path:

```bash
streamlit run app/streamlit_app.py
```
