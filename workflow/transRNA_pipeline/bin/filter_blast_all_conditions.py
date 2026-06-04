#!/usr/bin/env python3
"""
filter_blast_all_conditions.py
───────────────────────────────
Filters an annotated BLAST TSV (output of annotate_blast_lineage.py).

taxonkit reformat2 appends ONE semicolon-delimited lineage column as the
LAST column of each row:
  Domain;Kingdom;Phylum;Class;Order;Family;Genus;Species

Filtering criteria (ALL must pass):
  1.  Domain not in: Viruses, Virus, Acellular root, Acellular organisms
  2.  Phylum not in: excluded set from invertebrate_phyla.txt
  3.  Species contains no host fragment:
        homo sapiens, mus musculus, canis lupus, felis catus, apis mellifera
  4.  mismatch == 0
  5.  gapopen == 0
  6.  abs(qend - qstart) + 1 == alignment_length  (full-span)
  7.  alignment_length >= min_len  (default 18)

The lineage column is passed through to the output unchanged so the
taxonomy parser can re-use it directly.

Usage:
  python filter_blast_all_conditions.py [--min-len N] [--stats-out file]
      phyla.txt input_annotated.tsv output.tsv
"""

from pathlib import Path
import argparse
import re
import sys

COMMON_NAMES = {
    "sponges","coral","jellyfish","anemones","comb jellies","flatworms",
    "jaw worms","mesozoa","proboscis worms","gastrotrichs","rotifers",
    "roundworms","horsehair worms","mud dragons","spiny-crown worms",
    "acanthocephalans","spiny-headed worms","brush heads","pandora",
    "cycliophorans","goblet worms","marine mats","moss animals","bryozoans",
    "horseshoe worms","brachipods","lampshells","molluscs","slugs","snails",
    "squid","peanut worms","segmented worms","earthworms","ragworms",
    "spoon worms","beard worms","water bears","velvet worms","insects",
    "spiders","crabs","etc","starfish","urchins","arrow worms",
    "acorn worms","vertebrates","invertebrates",
}

VIRUS_DOMAINS = {"viruses","virus","acellular root","acellular organisms"}

HOST_FRAGS = [
    "homo sapiens","mus musculus","canis lupus","felis catus","apis mellifera"
]


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def clean_name(s):
    s = s.strip().strip(" ,;:.()[]{}")
    return re.sub(r"\s+", " ", s)


def extract_names(line):
    parts = re.split(r"[()]", line.strip())
    names = set()
    for part in parts:
        part = part.replace(",", " or ")
        for token in re.split(r"\s+or\s+", part):
            name = clean_name(token)
            if name and name.lower() not in COMMON_NAMES:
                names.add(name)
    return names


def load_excluded_phyla(path):
    names = set()
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                names.update(extract_names(line))
    return names


def is_host(species: str) -> bool:
    s = species.strip().lower()
    return any(f in s for f in HOST_FRAGS)


def parse_lineage(lineage_col: str) -> list:
    """
    Parse the semicolon-delimited lineage column (last column of annotated file).
    Returns list of 8 rank strings: [Domain, Kingdom, Phylum, ..., Species]
    """
    parts = [x.strip() for x in lineage_col.split(";")]
    while len(parts) < 8:
        parts.append("")
    return parts[:8]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("exclude_file")
    ap.add_argument("input_tsv")
    ap.add_argument("output_tsv")
    ap.add_argument("--min-len",   type=int, default=18)
    ap.add_argument("--stats-out", default=None)
    args = ap.parse_args()

    excluded = load_excluded_phyla(args.exclude_file)
    if not excluded:
        sys.exit("ERROR: No phyla names parsed from exclude file")

    total = kept = 0
    rm_v = rm_p = rm_h = rm_mm = rm_gap = rm_span = rm_short = bad = 0
    seen_queries = set()

    with open(args.input_tsv) as fin, open(args.output_tsv, "w") as fout:
        for line in fin:
            total += 1
            raw = line.rstrip("\n")
            fields = raw.split("\t")

            # Need at least 8 standard BLAST cols + 1 lineage col = 9 minimum
            if len(fields) < 9:
                bad += 1
                continue
            seen_queries.add(fields[0].strip())

            # Standard BLAST columns (0-based)
            try:
                aln_len  = int(fields[3])
                mismatch = int(fields[4])
                gapopen  = int(fields[5])
                qstart   = int(fields[6])
                qend     = int(fields[7])
            except ValueError:
                bad += 1
                continue

            qspan = abs(qend - qstart) + 1

            # Lineage: last column, semicolon-delimited
            ranks   = parse_lineage(fields[-1])
            domain  = ranks[0].lower()
            phylum  = ranks[2]
            species = ranks[7]

            # ── Filters ───────────────────────────────────────────────────────
            if domain in VIRUS_DOMAINS:    rm_v    += 1; continue
            if phylum in excluded:         rm_p    += 1; continue
            if is_host(species):           rm_h    += 1; continue
            if mismatch != 0:              rm_mm   += 1; continue
            if gapopen  != 0:              rm_gap  += 1; continue
            if qspan != aln_len:           rm_span += 1; continue
            if aln_len < args.min_len:     rm_short+= 1; continue

            fout.write(raw + "\n")
            kept += 1

    summary = (
        f"{Path(args.input_tsv).name}: total={total} kept={kept} | "
        f"rm_virus={rm_v} rm_phylum={rm_p} rm_host={rm_h} "
        f"rm_mismatch={rm_mm} rm_gap={rm_gap} "
        f"rm_not_fullspan={rm_span} rm_short={rm_short} "
        f"bad_format={bad} | excluded_phyla={len(excluded)}"
    )
    log(summary)

    # Write stats file
    stats_path = args.stats_out
    if stats_path is None:
        stats_path = str(Path(args.output_tsv).with_suffix("")) + "_filter_stats.tsv"
    stem = Path(args.input_tsv).name.replace("_blast_annotated.tsv", "")
    unique_queries = len(seen_queries)
    with open(stats_path, "w") as sf:
        sf.write("sample\ttool\ttotal\tunique_queries\tkept\t"
                 "rm_virus\trm_phylum\trm_host\trm_mismatch\t"
                 "rm_gap\trm_not_fullspan\trm_short\tbad\n")
        sf.write(
            f"{stem}\tblast\t{total}\t{unique_queries}\t{kept}\t"
            f"{rm_v}\t{rm_p}\t{rm_h}\t{rm_mm}\t"
            f"{rm_gap}\t{rm_span}\t{rm_short}\t{bad}\n"
        )


if __name__ == "__main__":
    main()