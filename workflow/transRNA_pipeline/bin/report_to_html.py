#!/usr/bin/env python3
"""
report_to_html.py
Converts pipeline_read_counts_report.tsv to styled HTML.
"""
import sys
import pandas as pd
from pathlib import Path

CSS = """
body{font-family:Arial,sans-serif;margin:24px;background:#f8f9fa;font-size:12px}
h1{color:#2c3e50;font-size:18px}
h2{color:#34495e;margin-top:20px;font-size:14px}
p.note{font-size:11px;color:#888;margin-top:4px}
table{border-collapse:collapse;margin-bottom:20px;width:100%}
th{background:#2c3e50;color:#fff;padding:6px 10px;text-align:right;white-space:nowrap;font-size:11px}
th:first-child{text-align:left}
td{padding:4px 10px;border-bottom:1px solid #dde;text-align:right;white-space:nowrap}
td:first-child{text-align:left;font-weight:bold}
tr:nth-child(even){background:#ecf0f1}
tr:hover{background:#d5dbdb}
.rj{color:#0070b8}
.st{color:#cc0000}
"""

COL_ORDER = [
    "Sample",
    "Raw reads",
    "Trimmed",
    "STAR mapped",
    "STAR unmapped",
    "Decon-a: Host+Human excl.",
    "Decon-b: 12 viral excl.",
    "MAGs excl.",
    "Candidate transRNAs total reads",
    "Candidate transRNAs unique sequences",
    "Candidate transRNAs ge5 duplicates",
    "Kraken classified total",
    "BLAST classified total",
    "Classified total",
    "Total transRNAs after all filters",
]

def fmt(x):
    try:
        v = int(x)
        return f"{v:,}" if v >= 0 else ""
    except Exception:
        return "" if (str(x) in ("nan","")) else str(x)

def group_of(sample):
    s = sample.upper()
    if s.startswith("RJ"): return "RJ"
    if s.startswith("T"):  return "ST"
    return ""

def make_table(df):
    present = [c for c in COL_ORDER if c in df.columns]
    th = "".join(f"<th>{c}</th>" for c in present)
    rows_html = []
    for _, row in df.iterrows():
        sample = str(row.get("Sample",""))
        grp = group_of(sample)
        cls = "rj" if grp=="RJ" else "st" if grp=="ST" else ""
        cells = []
        for col in present:
            val = fmt(row.get(col,""))
            td_cls = f' class="{cls}"' if col=="Sample" else ""
            cells.append(f"<td{td_cls}>{val}</td>")
        rows_html.append("<tr>" + "".join(cells) + "</tr>")
    return (
        f"<table><thead><tr>{th}</tr></thead>"
        f"<tbody>{''.join(rows_html)}</tbody></table>"
    )

def main():
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} report.tsv report.html")
    df = pd.read_csv(sys.argv[1], sep="\t")
    # sort RJ first, then ST
    def sort_key(s):
        s = str(s).upper()
        return (0 if s.startswith("RJ") else 1, s)
    df = df.iloc[df["Sample"].map(sort_key).argsort()]

    table_html = make_table(df)
    html = f"""<!DOCTYPE html>
<html><head><meta charset="UTF-8">
<title>transRNA Pipeline — Read Count Report</title>
<style>{CSS}</style>
</head><body>
<h1>transRNA Pipeline — Read Count Report</h1>
<p class="note">
  Decon-a/b/MAGs columns show reads REMAINING after each decontamination step.<br>
  Candidate totals use re-inflated counts (duplicate weights applied).<br>
  STAR mapped = Trimmed − STAR unmapped.
</p>
{table_html}
</body></html>"""
    Path(sys.argv[2]).write_text(html, encoding="utf-8")
    print(f"HTML written: {sys.argv[2]}", file=sys.stderr)

if __name__ == "__main__":
    main()