import os
from primer_designer.export_manager import export_csv, export_fasta, export_json, export_txt, export_pdf, ExportOptions


def test_export_json_txt_csv_fasta(tmp_path):
    options = ExportOptions(output_dir=str(tmp_path), prefix="test_report")
    data = {
        "summary": {"selected_forward": "ATGCGT...", "selected_reverse": "CGTACG..."},
        "scores": [100, 95, 90],
    }
    primers = [
        {"id": "FWD-001", "sequence": "ATGCGTACGTTAGC"},
        {"id": "REV-001", "sequence": "CGTACGTTAGCTAG"},
    ]
    json_path = export_json(data, options)
    txt_path = export_txt(data, options)
    csv_path = export_csv(primers, options)
    fasta_path = export_fasta(primers, options)

    assert os.path.exists(json_path)
    assert os.path.exists(txt_path)
    assert os.path.exists(csv_path)
    assert os.path.exists(fasta_path)

    with open(json_path, "r", encoding="utf-8") as fh:
        assert "selected_forward" in fh.read()
    with open(txt_path, "r", encoding="utf-8") as fh:
        assert "selected_forward" in fh.read()
    with open(csv_path, "r", encoding="utf-8") as fh:
        assert "FWD-001" in fh.read()
    with open(fasta_path, "r", encoding="utf-8") as fh:
        assert ">FWD-001" in fh.read()


def test_export_pdf_optional(tmp_path):
    options = ExportOptions(output_dir=str(tmp_path), prefix="test_pdf", include_pdf=True)
    data = {"summary": {"primer": "ATGCGT"}}
    pdf_path = export_pdf(data, options)
    if pdf_path is not None:
        assert os.path.exists(pdf_path)
        assert pdf_path.endswith(".pdf")
    else:
        assert pdf_path is None
