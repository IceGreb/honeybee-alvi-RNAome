#!/usr/bin/env python3
"""
collapse_stats.py
──────────────────
Reads the merged duplicated.detail.txt and writes a collapse_stats.tsv summary.

Usage:
  python3 collapse_stats.py <sample> <min_occ> <merged_reads>

  merged_reads: total reads in merged.fq before rmdup (passed from COLLAPSE_READS)
"""
import sys
from pathlib import Path


def count_groups(path: Path, min_occ: int):
    total = ge5 = 0
    with path.open() as f:
        for line in f:
            parts = line.strip().split("\t", 1)
            if not parts:
                continue
            try:
                c = int(parts[0].strip())
                total += 1
                if c >= min_occ:
                    ge5 += 1
            except ValueError:
                pass
    return total, ge5


def main():
    if len(sys.argv) < 3:
        sys.exit(f"Usage: {sys.argv[0]} <sample> <min_occ> [merged_reads]")

    sample       = sys.argv[1]
    min_occ      = int(sys.argv[2])
    merged_reads = int(sys.argv[3]) if len(sys.argv) > 3 else 0

    merged_path = Path(f"{sample}_merged_duplicated.detail.txt")

    if not merged_path.exists():
        print(f"Warning: no duplicated.detail.txt found for {sample}", file=sys.stderr)
        total, ge5 = 0, 0
    else:
        total, ge5 = count_groups(merged_path, min_occ)
        print(f"  {sample}: merged_reads={merged_reads}, "
              f"{total} total groups, {ge5} with >= {min_occ} duplicates",
              file=sys.stderr)

    out_path = Path(f"{sample}_collapse_stats.tsv")
    with out_path.open("w") as out:
        out.write("sample\tmate\tmerged_reads\ttotal_groups\tge5_groups\n")
        for mate in ["1", "2"]:
            out.write(f"{sample}\t{mate}\t{merged_reads}\t{total}\t{ge5}\n")

    print(f"Written: {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()