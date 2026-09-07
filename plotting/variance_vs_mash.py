
"""
Parses a binned fragments file where pairs of FASTA paths are followed by a
shared line of bin numbers. Each pair is treated as one plasmid pair data point.
Looks up the mash distance for each pair from a mash distance file (tab-separated,
distance in column 3), matching by filename after replacing 'plasmid_references'
with 'binned_fragments'. Produces a scatter plot of mash distance (x) vs bin
number variance (y), one point per plasmid pair, with Spearman rank correlation
and p-value.
 
Expected binned file structure:
    /binned_fragments/CAV1321_all_references_6_binned.fasta/binned_fragments/CAV1344_all_references_3_binned.fasta
    6 6 4 7 7 6 6 6 ...
 
Expected mash distance file structure (tab-separated):
    ./plasmid_references/A.fasta    ./plasmid_references/B.fasta    0.012    0    950/1000
"""
 
import sys
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
 
 
# ── helpers ───────────────────────────────────────────────────────────────────
 
def parse_binned_file(filepath):
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
 
 
def normalise_path(path):
    """
    Strip leading ./ and directory, return just the basename with
    'plasmid_references' replaced by 'binned_fragments' for matching.
    """
    path = path.strip()
    basename = os.path.basename(path)
    return basename
 
 
def parse_mash_file(filepath):
    """
    Parse mash distance file into a dict:
        (basename_a, basename_b) -> mash_distance
    Stored in both orders for easy lookup.
    'plasmid_references' is replaced with 'binned_fragments' in basenames.
    """
    distances = {}
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                continue
            try:
                dist = float(cols[2])
            except ValueError:
                continue
            a = normalise_path(cols[0]).replace("plasmid_references", "binned_fragments")
            b = normalise_path(cols[1]).replace("plasmid_references", "binned_fragments")
            distances[(a, b)] = dist
            distances[(b, a)] = dist  # store both orders
    
    return distances
 
 
def pair_label(paths):
    ids = []
    for p in paths:
        m = re.search(r'(CAV\d+)', os.path.basename(p))
        if m:
            ids.append(m.group(1))
    return " / ".join(ids) if ids else os.path.basename(paths[0])
 
 
# ── plotting ──────────────────────────────────────────────────────────────────
 
def plot_scatter(points):
    x = np.array([p["mash_distance"] for p in points])
    y = np.array([p["variance"] for p in points])
    labels = [p["label"] for p in points]
 
    fig, ax = plt.subplots(figsize=(9, 6))
 
    scatter = ax.scatter(
        x, y,
        s=90, zorder=3,
        c=y, cmap="viridis",
        edgecolors="white", linewidths=0.6
    )
 
    # for label, xi, yi in zip(labels, x, y):
    #     ax.annotate(
    #         label, (xi, yi),
    #         textcoords="offset points", xytext=(8, 4),
    #         fontsize=8
    #     )
 
    cbar = fig.colorbar(scatter, ax=ax, pad=0.02)
    cbar.set_label("Variance", fontsize=9)
 
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
        #     fontsize=9,
        #     verticalalignment="top",
        #     bbox=dict(boxstyle="round,pad=0.4", facecolor="lightyellow", edgecolor="gray", alpha=0.8)
        # )
        print(f"Spearman rank correlation (ρ): {r:.4f}")
        print(f"p-value: {p_value:.4f} ({sig_label} at α=0.05)")
        print(f"Note: based on subsample of {len(x)} pairs — interpret p-value with caution.")
    elif len(x) == 2:
        r, _ = stats.spearmanr(x, y)
        # ax.text(
        #     0.05, 0.95, f"Spearman ρ = {float(r):.4f}\n(need ≥3 pairs for p-value)",
        #     transform=ax.transAxes,
        #     fontsize=9, color="#ffa657",
        #     verticalalignment="top",
        #     bbox=dict(boxstyle="round,pad=0.4", facecolor="#1a1d27", edgecolor="#ffa657", alpha=0.8)
        # )
        print(f"Spearman ρ = {float(r):.4f} (need at least 3 pairs to compute p-value)")
    else:
        print("Not enough points to compute correlation (need at least 2).")
 
    ax.set_xlabel("Mash distance between plasmid pair", fontsize=11)
    ax.set_ylabel("Variance of distance distribution", fontsize=11)
   # ax.set_title("Plasmid pair: mash distance vs bin variance", fontsize=13)
    ax.grid(True, color="#30363d", linewidth=0.5, linestyle="--")
 
    plt.tight_layout()
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "binned_analysis.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Plot saved to: {out_path}")
    plt.show()
 
 
# ── main ──────────────────────────────────────────────────────────────────────
 
def main():
    if len(sys.argv) != 3:
        print("Usage: python analyze_binned_fragments.py <binned_file> <mash_distance_file>")
        sys.exit(1)
 
    binned_file = sys.argv[1]
    mash_file = sys.argv[2]
 
    for path, label in [(binned_file, "binned file"), (mash_file, "mash distance file")]:
        if not os.path.isfile(path):
            print(f"Error: {label} not found: {path}")
            sys.exit(1)
 
    print(f"Parsing binned file: {binned_file}")
    entries = parse_binned_file(binned_file)
    print(f"Found {len(entries)} plasmid pair(s).\n")
 
    print(f"Parsing mash distance file: {mash_file}")
    mash_distances = parse_mash_file(mash_file)
    print(f"Loaded {len(mash_distances) // 2} unique mash distance entries.\n")
 
    points = []
    for fasta_paths, numbers in entries:
        #label = pair_label(fasta_paths)
        print(label)
        if len(fasta_paths) < 2:
            print(f"  WARNING: pair [{label}] has fewer than 2 files. Skipping.")
            continue
 
        a = os.path.basename(fasta_paths[0])
        b = os.path.basename(fasta_paths[1])

        a = a.replace("_binned","")
        b = b.replace("_binned","")
        dist = mash_distances.get((a, b))
        if dist is None:
            print(f"  WARNING: no mash distance found for pair [{label}] ({a} vs {b}). Skipping.")
            continue
 
        variance = float(np.var(numbers))
        print(f"  Pair [{label}]: mash distance = {dist:.6f}, variance = {variance:.4f}")
        points.append({"label": label, "mash_distance": dist, "variance": variance})
 
    if not points:
        print("No valid plasmid pairs found. Cannot plot.")
        sys.exit(1)
 
    print(f"\nPlotting {len(points)} pairs...")
    plot_scatter(points)
 
 
if __name__ == "__main__":
    main()