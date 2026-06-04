#!/usr/bin/env python3
"""
make_reads_posttrim_tab.py
──────────────────────────
Reads all *_trimmed_stats.tsv files in the current directory (staged by
Nextflow) and writes reads_posttrim_tab.tsv with columns Sample, TotalReads.
TotalReads = mate-1 read count from seqkit stats, used as RPM denominator.

Usage (called from Nextflow process, no arguments needed):
  python3 make_reads_posttrim_tab.py
"""
import glob
import re
import sys
from pathlib import Path

import pandas as pd

rows = []
for f in sorted(glob.glob("*_trimmed_stats.tsv")):
    sample = re.sub(r"_trimmed_stats\.tsv$", "", Path(f).name)
    try:
        df = pd.read_csv(f, sep="\t")
    except Exception as e:
        print(f"Warning: could not read {f}: {e}", file=sys.stderr)
        continue

    # Prefer mate-1 row; fall back to first row if not identifiable
    mate1 = df[df["file"].astype(str).str.contains(r"val_1|_1\.fq|clean1|mate1", na=False, regex=True)]
    if mate1.empty:
        mate1 = df.head(1)

    if mate1.empty:
        print(f"Warning: no rows found in {f}", file=sys.stderr)
        continue

    total = int(mate1["num_seqs"].iloc[0])
    rows.append({"Sample": sample, "TotalReads": total})
    print(f"  {sample}: {total:,} reads (mate 1)", file=sys.stderr)

out = pd.DataFrame(rows, columns=["Sample", "TotalReads"])
out.to_csv("reads_posttrim_tab.tsv", sep="\t", index=False)
print(f"reads_posttrim_tab.tsv written: {len(out)} samples", file=sys.stderr)
