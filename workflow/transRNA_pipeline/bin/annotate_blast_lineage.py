#!/usr/bin/env python3
"""
annotate_blast_lineage.py
─────────────────────────
Annotates a BLAST tabular output file with an 8-rank lineage column.

Writes the temp taxid file to the current working directory (the Nextflow
work dir on RDS) instead of /tmp, to avoid filling the login node's /tmp.
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

TAXID_SEP = re.compile(r"[;,]")

LINEAGE_FMT = (
    "{domain|acellular root|superkingdom}"
    ";{kingdom};{phylum};{class};{order};{family};{genus};{species}"
)

FALLBACK = ";;;;;;;"   # 7 semicolons = 8 empty fields


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def first_taxid(field: str) -> int:
    for part in TAXID_SEP.split(field.strip()):
        part = part.strip()
        if part.isdigit():
            v = int(part)
            if v > 0:
                return v
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input",     required=True)
    ap.add_argument("--output",    required=True)
    ap.add_argument("--taxid-col", type=int, default=13)
    args = ap.parse_args()

    taxid_idx = args.taxid_col - 1

    log(f"Pass 1: reading {args.input} (taxid col={args.taxid_col})")
    lines = []
    line_taxids = []

    with open(args.input) as fh:
        for line in fh:
            raw = line.rstrip("\n")
            lines.append(raw)
            fields = raw.split("\t")
            tid = first_taxid(fields[taxid_idx]) if len(fields) > taxid_idx else 0
            line_taxids.append(tid)

    unique_taxids = sorted({t for t in line_taxids if t != 0})
    log(f"  {len(lines)} lines, {len(unique_taxids)} unique non-zero taxids")

    # ── Write temp taxid file to current dir (Nextflow work dir on RDS) ───────
    # Avoids filling /tmp on the login/compute node
    tmp_path = Path(".") / f"_taxids_tmp_{Path(args.output).stem}.txt"
    log(f"Running taxonkit reformat2 (temp file: {tmp_path})...")

    try:
        with tmp_path.open("w") as tmp:
            for tid in unique_taxids:
                tmp.write(f"{tid}\n")

        proc = subprocess.run(
            ["taxonkit", "reformat2", "-I", "1", "-f", LINEAGE_FMT, str(tmp_path)],
            capture_output=True,
            text=True,
        )
    finally:
        tmp_path.unlink(missing_ok=True)

    if proc.returncode != 0:
        log(f"[ERROR] taxonkit reformat2 failed (exit {proc.returncode}):")
        log(proc.stderr[:2000])
        sys.exit(1)

    if proc.stderr.strip():
        log(proc.stderr.strip())

    # ── Parse output: {taxid: lineage_string} ────────────────────────────────
    lineage_map = {}
    for out_line in proc.stdout.splitlines():
        parts = out_line.strip().split("\t")
        if len(parts) < 2:
            continue
        try:
            tid = int(parts[0])
        except ValueError:
            continue
        lineage_map[tid] = parts[1].strip()

    log(f"  {len(lineage_map)} lineages resolved")

    # ── Pass 2: write annotated output ───────────────────────────────────────
    log(f"Pass 2: writing {args.output}")
    no_lineage = 0

    with open(args.output, "w") as out:
        for raw, tid in zip(lines, line_taxids):
            lineage = lineage_map.get(tid, FALLBACK)
            if not lineage:
                lineage = FALLBACK
                no_lineage += 1
            out.write(raw + "\t" + lineage + "\n")

    if no_lineage:
        log(f"  {no_lineage} lines used fallback (taxid 0 or not resolved)")
    log("Done.")


if __name__ == "__main__":
    main()