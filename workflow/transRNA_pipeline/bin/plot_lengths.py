#!/usr/bin/env python3
"""
plot_lengths.py
───────────────
Reads *_1_ge5_detected.hist files (weighted length histograms produced by the
ge5 script) and plots the average reinflated length distribution for RJ and ST
sample groups over the 15–200 nt window.

Hist file format (no header):
    length<TAB>weighted_count

Only mate-1 hist files are used (*_1_ge5_detected.hist); both mates carry
identical content in the merged-pool pipeline, so mate 2 is skipped to avoid
double-counting.

Sample grouping:
    R* → RJ
    T* → ST

Output (one plot):
    <date>_avg_reinflated_length_distribution_15to100nt.png
"""

import matplotlib
matplotlib.use("Agg")

import sys
from datetime import date
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PLOT_MIN = 15
PLOT_MAX = 100
HIGHLIGHT_TICKS = [20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100]


def load_hist(path: Path) -> pd.DataFrame:
    """Load a length-histogram file and return a percent-normalised DataFrame."""
    df = pd.read_csv(path, sep="\t", header=None, names=["length", "count"])
    df = df.dropna()
    df["length"] = df["length"].astype(int)
    df["count"]  = pd.to_numeric(df["count"], errors="coerce").fillna(0)
    total = df["count"].sum()
    if total <= 0:
        return pd.DataFrame(columns=["length", "percent"])
    df["percent"] = df["count"] / total * 100.0
    return df[["length", "percent"]]


def sample_name(path: Path) -> str:
    """Strip '_ge5_detected.hist' suffix to get the base sample name."""
    return path.name.replace("_ge5_detected.hist", "")


def infer_group(name: str) -> str:
    if name.upper().startswith("R"):
        return "RJ"
    if name.upper().startswith("T"):
        return "ST"
    return None


def align_and_mean(dfs: list, bins: np.ndarray) -> np.ndarray:
    """Align each sample to a shared bin axis and return the mean curve."""
    if not dfs:
        return np.zeros(len(bins))
    aligned = []
    for df in dfs:
        # Restrict to window and re-normalise within it
        sub = df[df["length"].between(PLOT_MIN, PLOT_MAX)].copy()
        sub_total = sub["percent"].sum()
        if sub_total > 0:
            sub["percent"] = sub["percent"] / sub_total * 100.0
        s = sub.set_index("length")["percent"].reindex(bins, fill_value=0.0)
        aligned.append(s.to_numpy())
    return np.vstack(aligned).mean(axis=0)


def make_plot(means: dict, bins: np.ndarray, grouped: dict, outfile: Path):
    fig, ax = plt.subplots(figsize=(12, 7))

    # Alternating column shading
    for i, b in enumerate(bins):
        if i % 2 == 1:
            ax.axvspan(b - 0.5, b + 0.5, color="gray", alpha=0.08, zorder=0)

    # Vertical dotted guide lines
    for x in HIGHLIGHT_TICKS:
        if PLOT_MIN <= x <= PLOT_MAX:
            ax.axvline(x=x, color="gray", linestyle=":", linewidth=1.8,
                       alpha=1.0, zorder=1)

    if grouped.get("ST"):
        ax.plot(bins, means["ST"], color="#ff0000", lw=2, label="ST", zorder=3)
    if grouped.get("RJ"):
        ax.plot(bins, means["RJ"], color="#00a4ff", lw=2, label="RJ", zorder=3)

    ax.set_xlabel("Sequence length (nt)", fontsize=13)
    ax.set_ylabel("Weighted percentage of reads (%)", fontsize=13)
    ax.set_title(f"Average Length Distribution ({PLOT_MIN}–{PLOT_MAX} nt)",
                 fontsize=15)
    ax.legend(frameon=True, facecolor="white")
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    ax.set_ylim(bottom=0)
    ax.set_xlim(PLOT_MIN - 0.5, PLOT_MAX + 0.5)

    ticks = sorted({PLOT_MIN} | {t for t in HIGHLIGHT_TICKS if PLOT_MIN <= t <= PLOT_MAX})
    ax.set_xticks(ticks)
    ax.tick_params(axis="x", rotation=45)

    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()
    print(f"Saved: {outfile}", file=sys.stderr)


def main():
    hist_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
    today    = date.today().strftime("%d_%m_%Y")

    hist_files = sorted(hist_dir.glob("*_ge5_detected.hist"))
    if not hist_files:
        print(f"No *_1_ge5_detected.hist files found in {hist_dir}", file=sys.stderr)
        sys.exit(1)

    grouped = {"RJ": [], "ST": []}

    for f in hist_files:
        name  = sample_name(f)
        group = infer_group(name)
        if group is None:
            print(f"Skipping {f.name}: cannot infer group from name '{name}'",
                  file=sys.stderr)
            continue

        df = load_hist(f)
        if df.empty or df["percent"].sum() == 0:
            print(f"Skipping {f.name}: empty or zero-weight histogram",
                  file=sys.stderr)
            continue

        grouped[group].append(df)
        print(f"Loaded {f.name} → group={group}, "
              f"n_bins={len(df)}, sum%={df['percent'].sum():.4f}",
              file=sys.stderr)

    if not grouped["RJ"] and not grouped["ST"]:
        print("No valid hist files found for RJ or ST.", file=sys.stderr)
        sys.exit(1)

    bins  = np.arange(PLOT_MIN, PLOT_MAX + 1)
    means = {label: align_and_mean(grouped[label], bins)
             for label in ["RJ", "ST"]}

    for label in ["RJ", "ST"]:
        if grouped[label]:
            print(f"{label}: n={len(grouped[label])}, "
                  f"mean curve sum={means[label].sum():.4f}%",
                  file=sys.stderr)

    outfile = hist_dir / f"{today}_avg_reinflated_length_distribution_15to100nt.png"
    make_plot(means, bins, grouped, outfile)


if __name__ == "__main__":
    main()
