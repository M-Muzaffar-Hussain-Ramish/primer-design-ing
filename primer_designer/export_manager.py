import csv
import json
import os
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

@dataclass
class ExportOptions:
    output_dir: str = "results"
    prefix: str = "report"
    include_pdf: bool = False


def _ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def export_json(data: Dict[str, Any], options: ExportOptions) -> str:
    _ensure_dir(options.output_dir)
    path = os.path.join(options.output_dir, f"{options.prefix}.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2)
    return path


def export_txt(data: Dict[str, Any], options: ExportOptions) -> str:
    _ensure_dir(options.output_dir)
    path = os.path.join(options.output_dir, f"{options.prefix}.txt")
    lines: List[str] = []
    for key, value in data.items():
        lines.append(f"{key}:")
        if isinstance(value, dict):
            for k, v in value.items():
                lines.append(f"  {k}: {v}")
        elif isinstance(value, list):
            for item in value:
                lines.append(f"  - {item}")
        else:
            lines.append(f"  {value}")
        lines.append("")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    return path


def export_csv(primer_table: List[Dict[str, Any]], options: ExportOptions) -> str:
    _ensure_dir(options.output_dir)
    path = os.path.join(options.output_dir, f"{options.prefix}.csv")
    if not primer_table:
        with open(path, "w", encoding="utf-8", newline="") as fh:
            fh.write("")
        return path
    fieldnames = list(primer_table[0].keys())
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in primer_table:
            writer.writerow({k: row.get(k, "") for k in fieldnames})
    return path


def export_fasta(primers: List[Dict[str, Any]], options: ExportOptions) -> str:
    _ensure_dir(options.output_dir)
    path = os.path.join(options.output_dir, f"{options.prefix}.fasta")
    with open(path, "w", encoding="utf-8") as fh:
        for p in primers:
            seq = p.get("sequence", "")
            header = p.get("id", "primer")
            fh.write(f">{header}\n")
            fh.write(seq + "\n")
    return path


def export_pdf(data: Dict[str, Any], options: ExportOptions) -> Optional[str]:
    if not options.include_pdf:
        return None
    try:
        from reportlab.lib.pagesizes import letter
        from reportlab.pdfgen import canvas
    except ImportError:
        return None

    _ensure_dir(options.output_dir)
    path = os.path.join(options.output_dir, f"{options.prefix}.pdf")
    c = canvas.Canvas(path, pagesize=letter)
    width, height = letter
    y = height - 40
    c.setFont("Helvetica-Bold", 14)
    c.drawString(40, y, f"Primer Design Report: {options.prefix}")
    y -= 24
    c.setFont("Helvetica", 10)
    for key, value in data.items():
        if y < 40:
            c.showPage()
            y = height - 40
            c.setFont("Helvetica", 10)
        c.drawString(40, y, f"{key}:")
        y -= 14
        if isinstance(value, dict):
            for k, v in value.items():
                if y < 40:
                    c.showPage()
                    y = height - 40
                    c.setFont("Helvetica", 10)
                c.drawString(60, y, f"{k}: {v}")
                y -= 12
        elif isinstance(value, list):
            for item in value:
                if y < 40:
                    c.showPage()
                    y = height - 40
                    c.setFont("Helvetica", 10)
                c.drawString(60, y, f"- {item}")
                y -= 12
        else:
            c.drawString(60, y, str(value))
            y -= 12
        y -= 10
    c.save()
    return path
