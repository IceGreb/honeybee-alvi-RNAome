#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import argparse
import gzip
import sys


def log(msg):
    print(msg, file=sys.stderr)


def open_maybe_gz(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def iter_fastq(path):
    with open_maybe_gz(path) as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline().rstrip("\n")
            p = fh.readline()
            q = fh.readline()
            if not q:
                raise ValueError(f"Incomplete FASTQ record in {path}")
            yield h.rstrip("\n"), s


def read_seq_lengths(fq_path):
    """Return {seqid: length} from a FASTQ file (representative collapsed seqs)."""
    lengths = {}
    for header, seq in iter_fastq(fq_path):
        seqid = header.lstrip("@").split()[0]
        lengths[seqid] = len(seq)
    return lengths


def read_dup_weights(path, min_occ=5):
    """
    Parse duplicated.detail.txt (format: count<TAB>ID1,ID2,...).
    Returns {first_id: count} for groups with count >= min_occ.
    Only the first ID per group is used — it is the representative sequence
    present in the collapsed-clean FASTQ.
    """
    weights = {}
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t", 1)
            if len(parts) != 2:
                continue
            try:
                occ = int(parts[0].strip())
            except ValueError:
                continue
            if occ < min_occ:
                continue
            first_id = parts[1].split(",")[0].strip()
            if first_id:
                weights[first_id] = occ
    return weights


def write_ids(path, ids):
    with path.open("w") as out:
        for seqid in sorted(ids):
            out.write(f"{seqid}\n")


def write_weighted_ids(path, kept):
    with path.open("w") as out:
        for seqid in sorted(kept):
            length, weight, source = kept[seqid]
            out.write(f"{seqid}\t{length}\t{weight}\t{source}\n")


def write_weighted_hist(path, length_weight_pairs):
    hist = Counter()
    for length, weight in length_weight_pairs:
        hist[length] += weight
    with path.open("w") as out:
        for length in sorted(hist):
            out.write(f"{length}\t{hist[length]}\n")


def run_filtered_mode(args):
    dup_file = Path(args.dup_file)
    fq_file  = Path(args.fq_file)
    sample   = args.sample
    out_dir  = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Lengths from collapsed-clean FASTQ (representatives only)
    seq_lengths = read_seq_lengths(fq_file)

    # Sequences passing >= min_occ threshold
    weights = read_dup_weights(dup_file, min_occ=args.min_occ)

    # Apply min_len filter
    kept = {}
    for seqid, weight in sorted(weights.items()):
        length = seq_lengths.get(seqid, 0)
        if length >= args.min_len:
            kept[seqid] = (length, weight, "merged")

    log(
        f"{sample}"
        f"\tinput_ids={len(weights)}"
        f"\tinput_weight={sum(weights.values())}"
        f"\tkept_ids={len(kept)}"
        f"\tkept_weight={sum(v[1] for v in kept.values())}"
    )

    write_ids(out_dir / f"{sample}_ge5_detected_ids.txt", kept.keys())
    write_weighted_ids(out_dir / f"{sample}_ge5_detected_weighted_ids.tsv", kept)
    write_weighted_hist(
        out_dir / f"{sample}_ge5_detected.hist",
        [(length, weight) for length, weight, source in kept.values()],
    )


def fastq_basename(path):
    name = path.name
    for suffix in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def run_fastq_mode(args):
    fastq_dir = Path(args.fastq_dir)
    out_dir   = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(fastq_dir.glob(args.fastq_pattern))
    if not files:
        sys.stderr.write(
            f"No FASTQ files found in {fastq_dir} with pattern {args.fastq_pattern}\n"
        )
        sys.exit(1)

    for fq in files:
        hist = Counter()
        n = 0
        for header, seq in iter_fastq(fq):
            qlen = len(seq)
            if qlen >= args.min_len:
                hist[qlen] += 1
                n += 1

        out_path = out_dir / f"{fastq_basename(fq)}.hist"
        with out_path.open("w") as out:
            for length in sorted(hist):
                out.write(f"{length}\t{hist[length]}\n")

        log(f"{fq.name}\treads_kept={n}\tlength_bins={len(hist)}\toutput={out_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["filtered", "fastq"], default="filtered")

    # filtered mode
    ap.add_argument("--dup-file",    help="Path to *_merged_duplicated.detail.txt")
    ap.add_argument("--fq-file",     help="Path to *_merged_collapsed_clean.fq")
    ap.add_argument("--sample",      help="Sample name (used for output filenames)")

    # fastq mode
    ap.add_argument("--fastq-dir")
    ap.add_argument("--fastq-pattern", default="*.fq*")

    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--min-occ",  type=int, default=5)
    ap.add_argument("--min-len",  type=int, default=18)

    args = ap.parse_args()

    if args.mode == "filtered":
        for flag, name in [("--dup-file", args.dup_file),
                           ("--fq-file",  args.fq_file),
                           ("--sample",   args.sample)]:
            if not name:
                sys.stderr.write(f"{flag} is required in --mode filtered\n")
                sys.exit(1)
        run_filtered_mode(args)

    elif args.mode == "fastq":
        if not args.fastq_dir:
            sys.stderr.write("--fastq-dir is required in --mode fastq\n")
            sys.exit(1)
        run_fastq_mode(args)


if __name__ == "__main__":
    main()
