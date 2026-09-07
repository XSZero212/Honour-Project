
"""
Scatterplot: distance (file 1) vs. third-column ratio (file 2)
 
File 1 format (TSV):  plasmid_1  plasmid_2  distance
File 2 format (TSV):  plasmid_1  plasmid_2  col3  col4  col5
  where col3 is a fraction like  1000/1000
# """
 
import argparse
import sys
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy import stats
 
# ── helpers ──────────────────────────────────────────────────────────────────
 
def parse_fraction(value: str) -> float:
    """Convert '950/1000' → 0.95.  Plain numbers pass through unchanged."""
    value = str(value).strip()
    if "/" in value:
        num, den = value.split("/", 1)
        return float(num) / float(den)
    return float(value)
 
 
def make_key(a: str, b: str) -> frozenset:
    """Order-independent pair key so AP1–AP2 matches AP2–AP1."""
    return frozenset([a.strip(), b.strip()])
 
 
# ── loaders ──────────────────────────────────────────────────────────────────
 
def load_distances(path: str) -> dict:
    """Returns {frozenset(p1, p2): distance}."""
    df = pd.read_csv(path, sep="\t", header=0,
                     names=["plasmid_1", "plasmid_2", "distance"])
    return {make_key(r.plasmid_1, r.plasmid_2): float(r.distance)
            for _, r in df.iterrows()}
 
 
def load_ratios(path: str) -> dict:
    """Returns {frozenset(p1, p2): ratio} using the 3rd data column."""
    df = pd.read_csv(path, sep="\t", header=None)
    # Strip .fna suffixes if present so names match the distance file
    df[0] = df[0].str.replace(r"\.fna$", "", regex=True)
    df[1] = df[1].str.replace(r"\.fna$", "", regex=True)
    print(df)
    return {make_key(r[0], r[1]): float(r[2])
            for _, r in df.iterrows() if r[2]<=0.15}
 
def parse_mash_file(filepath):
    """Returns dict: (basename_a, basename_b) -> mash_distance (both orders)."""
    distances = []
    with open(filepath) as fh:
        for line in fh:
            cols = line.strip().split("\t")
            if len(cols) < 3:
                continue
            try:
                dist = float(cols[2])
            except ValueError:
                continue
            a = os.path.basename(cols[0]).replace("plasmid_references", "binned_fragments")
            b = os.path.basename(cols[1]).replace("plasmid_references", "binned_fragments")
            distances.append = (frozenset(a,b),dist)
            distances.append = (frozenset(a,b),dist)
    return distances
# ── main ─────────────────────────────────────────────────────────────────────
 
def main():
    parser = argparse.ArgumentParser(
        description="Scatterplot: pairwise distance vs. alignment ratio"
    )
    parser.add_argument("distance_file",
                        help="TSV with columns: plasmid_1, plasmid_2, distance")
    parser.add_argument("ratio_file",
                        help="TSV with columns: plasmid_1, plasmid_2, col3 (fraction), …")
    parser.add_argument("-o", "--output", default="scatterplot.png",
                        help="Output image path (default: scatterplot.png)")
    parser.add_argument("--xlabel", default="Pairwise distance",
                        help="X-axis label")
    parser.add_argument("--ylabel", default="Alignment ratio (col 3)",
                        help="Y-axis label")
    parser.add_argument("--title", default="Distance vs. Alignment ratio",
                        help="Plot title")
    parser.add_argument("--alpha", type=float, default=0.7,
                        help="Point transparency (0–1, default 0.7)")
    parser.add_argument("--point-size", type=float, default=40,
                        help="Marker size (default 40)")
    args = parser.parse_args()
 
    # Load both files
    distances = load_distances(args.distance_file)
    ratios    = load_ratios(args.ratio_file)
 
    # Inner join on shared pairs
    common_keys = set(distances) & set(ratios)
    if not common_keys:
        sys.exit(
            "ERROR: No matching plasmid pairs found between the two files.\n"
            "Check that names match (the script strips .fna suffixes automatically)."
        )
 
    x = [distances[k] for k in common_keys]
    y = [ratios[k]    for k in common_keys]
 
    print(f"Matched {len(common_keys)} pairs  "
          f"({len(distances) - len(common_keys)} distance-only, "
          f"{len(ratios) - len(common_keys)} ratio-only — skipped)")
 
    # ── statistics ────────────────────────────────────────────────────────────
    # Pearson (linear correlation)
    pearson_r, pearson_p = stats.pearsonr(x, y)
 
    # Spearman (rank-based, no normality assumption)
    spearman_r, spearman_p = stats.spearmanr(x, y)
 
    # Shapiro-Wilk normality test on both variables (only reliable for n < 5000)
    sw_stat_x, sw_p_x = stats.shapiro(x) if len(x) <= 5000 else (float("nan"), float("nan"))
    sw_stat_y, sw_p_y = stats.shapiro(y) if len(y) <= 5000 else (float("nan"), float("nan"))
 
    def fmt_p(p):
        if p < 0.001:
            return "< 0.001"
        return f"= {p:.3f}"
 
    print("\n── Statistical results ──────────────────────────────")
    print(f"  n pairs            : {len(x)}")
    print(f"  Pearson  r         : {pearson_r:+.4f}   p {fmt_p(pearson_p)}")
    print(f"  Spearman ρ         : {spearman_r:+.4f}   p {fmt_p(spearman_p)}")
    print(f"  Shapiro-Wilk (x)   : W = {sw_stat_x:.4f},  p {fmt_p(sw_p_x)}")
    print(f"  Shapiro-Wilk (y)   : W = {sw_stat_y:.4f},  p {fmt_p(sw_p_y)}")
    print("─────────────────────────────────────────────────────\n")
 
    # ── plot ──────────────────────────────────────────────────────────────────
    #sns.set_theme(style="whitegrid", context="notebook")
 
    fig, ax = plt.subplots(figsize=(8, 6))
 
    ax.scatter(x, y,
               s=args.point_size,
               alpha=args.alpha,
               edgecolors="white",
               linewidths=0.4,
               #color=sns.color_palette("mako", 1)[0]
               )
 
    # Regression line with 95 % CI
    # sns.regplot(x=x, y=y, ax=ax, scatter=False,
    #             line_kws={"color": "crimson", "linewidth": 1.5, "linestyle": "--"},
    #             ci=95)
 
    # Annotation box with both coefficients
    annotation = (
        f"Spearman ρ = {spearman_r:+.3f}  (p {fmt_p(spearman_p)})\n"
        f"n = {len(x)}"
    )
    # ax.text(0.03, 0.97, annotation,
    #         transform=ax.transAxes,
    #         fontsize=9.5,
    #         verticalalignment="top",
    #         bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
    #                   edgecolor="lightgray", alpha=0.9))
 
    ax.set_xlabel("Mash Distance", fontsize=13)
    ax.set_ylabel("Pling Distance", fontsize=13)
    ax.set_title("", fontsize=15, fontweight="bold", pad=12)
 
    plt.tight_layout()
    plt.savefig(args.output, dpi=150)
    print(f"Plot saved → {args.output}")
 
 
if __name__ == "__main__":
    main()