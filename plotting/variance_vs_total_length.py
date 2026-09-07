
"""
Parses a binned fragments file where pairs of FASTA paths are followed by a
shared line of bin numbers. Each pair is treated as one plasmid pair data point.
Computes the total sequence length of each FASTA file found under
'fragmented_database_final', averages across the two files in each pair, then
produces a scatter plot of average total length (x) vs bin number variance (y),
one point per plasmid pair, with Spearman rank correlation and p-value.
 
Expected input file structure (two paths, then their shared numbers, repeating):
    /binned_fragments/CAV1321_all_references_6_binned.fasta/binned_fragments/CAV1344_all_references_3_binned.fasta
    6 6 4 7 7 6 6 6 ...
    /binned_fragments/CAV1400_all_references_2_binned.fasta/binned_fragments/CAV1401_all_references_5_binned.fasta
    5 5 6 7 6 5 7 6 ...
"""
 
import sys
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
 
 
# ── helpers ───────────────────────────────────────────────────────────────────
 
def parse_binned_file(filepath):
    """
    Parse the binned file into a list of ([fasta_paths], [numbers]) tuples.
    Each group of FASTA paths (possibly concatenated on one line) is paired
    with the number sequence on the following line.
    """
    entries = []
    current_paths = []
 
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
 
            if ".fasta" in line:
                parts = re.split(r'(?=/)', line)
                for part in parts:
                    part = part.strip()
                    if part.endswith(".fasta"):
                        current_paths.append(part)
            else:
                tokens = line.split()
                numbers = []
                for t in tokens:
                    try:
                        numbers.append(int(t))
                    except ValueError:
                        pass
 
                if numbers and current_paths:
                    entries.append((list(current_paths), numbers))
                    current_paths = []
 
    return entries
 
 
def find_fasta_in_db(fasta_path, db_dir):
    """Walk db_dir to find a file matching the basename of fasta_path."""
    basename = os.path.basename(fasta_path)
    basename = basename.replace("_binned","")
    for root, dirs, files in os.walk(db_dir):
        if basename in files:
            return os.path.join(root, basename)
    return None
 
 
def total_sequence_length(filepath):
    """Sum the lengths of all sequences in a multi-FASTA file (ignores header lines)."""
    total = 0
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith(">"):
                total += len(line)
    return total
 
 
def pair_label(paths):
    """Build a short readable label from two FASTA basenames."""
    ids = []
    for p in paths:
        m = re.search(r'(CAV\d+)', os.path.basename(p))
        if m:
            ids.append(m.group(1))
    return " / ".join(ids) if ids else os.path.basename(paths[0])
 
 
# ── plotting ──────────────────────────────────────────────────────────────────
 
def plot_scatter(points):
    """
    points: list of dicts with keys 'label', 'avg_total_length', 'variance'
    """
    x = np.array([p["avg_total_length"] for p in points])
    y = np.array([p["variance"] for p in points])
    labels = [p["label"] for p in points]

    fig, ax = plt.subplots(figsize=(9, 6))

    scatter = ax.scatter(
        x, y,
        s=90, zorder=3,
        edgecolors="black", linewidths=0.6
    )

    # for label, xi, yi in zip(labels, x, y):
    #     ax.annotate(
    #         label, (xi, yi),
    #         textcoords="offset points", xytext=(8, 4),
    #         fontsize=8, color="black"
    #     )


    # Spearman rank correlation + significance test
    if len(x) >= 3:
        r, p_value = stats.spearmanr(x, y)
        r, p_value = float(r), float(p_value)
        sig_label = "significant" if p_value < 0.05 else "not significant"
        annotation = (
            f"Spearman ρ = {r:.4f}\n"
            f"p-value = {p_value:.4f}\n"
            f"({sig_label} at α=0.05)\n"
            f"n = {len(x)} pairs (subsample)"
        )
        # ax.text(
        #     0.05, 0.95, annotation,
        #     transform=ax.transAxes,
        #     fontsize=9, color="black",
        #     verticalalignment="top",
        #     bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="black", alpha=0.8)
        # )
        print(f"Spearman rank correlation (ρ): {r:.4f}")
        print(f"p-value: {p_value:.4f} ({sig_label} at α=0.05)")
        print(f"Note: based on subsample of {len(x)} pairs — interpret p-value with caution.")
    elif len(x) == 2:
        r, _ = stats.spearmanr(x, y)
        # ax.text(0.05, 0.95, annotation, transform=ax.transAxes, fontsize=9,
        #     verticalalignment="top",
        #     bbox=dict(boxstyle="round,pad=0.4", facecolor="lightyellow", edgecolor="gray", alpha=0.8))
        print(f"Spearman ρ = {float(r):.4f} (need at least 3 pairs to compute p-value)")
    else:
        print("Not enough points to compute correlation (need at least 2).")

    ax.set_xlabel("Average total plasmid length (bp)", fontsize=11)
    ax.set_ylabel("Variance of Pling distances", fontsize=11)
    #ax.set_title("Plasmid pair: avg total length vs bin variance", fontsize=13)
    ax.grid(True, color="lightgrey", linewidth=0.5, linestyle="--")

    plt.tight_layout()
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "binned_analysis.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Plot saved to: {out_path}")
    plt.show()
 
 
# ── main ──────────────────────────────────────────────────────────────────────
 
def main():
    if len(sys.argv) != 3:
        print("Usage: python analyze_binned_fragments.py <binned_file> <fragmented_database_final_dir>")
        sys.exit(1)
 
    binned_file = sys.argv[1]
    db_dir = sys.argv[2]
 
    if not os.path.isfile(binned_file):
        print(f"Error: binned file not found: {binned_file}")
        sys.exit(1)
    if not os.path.isdir(db_dir):
        print(f"Error: database directory not found: {db_dir}")
        sys.exit(1)
 
    print(f"Parsing: {binned_file}")
    entries = parse_binned_file(binned_file)
    print(f"Found {len(entries)} plasmid pair(s).")
    
    points = []
    for fasta_paths, numbers in entries:
        label = pair_label(fasta_paths)
        total_lengths = []
 
        for fp in fasta_paths:
            name = os.path.basename(fp)
            full_path = find_fasta_in_db(fp, db_dir)
            if full_path is None:
                print(f"  WARNING: '{name}' not found in '{db_dir}'. Skipping pair.")
                break
            length = total_sequence_length(full_path)
            print(f"  {name}: total length = {length:,} bp")
            total_lengths.append(length)
        else:
            avg_length = float(np.mean(total_lengths))
            variance = float(np.var(numbers))
            print(f"  Pair [{label}]: avg total length = {avg_length:,.1f} bp, variance = {variance:.4f}\n")
            points.append({"label": label, "avg_total_length": avg_length, "variance": variance})
 
    if not points:
        print("No valid plasmid pairs found. Cannot plot.")
        sys.exit(1)
 
    plot_scatter(points)
 
 
if __name__ == "__main__":
    main()