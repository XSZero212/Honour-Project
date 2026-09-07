"""
plot_vs_ground_trhuth.py — binned (real short-read) distances vs ground truth.

Compares the distance distribution sampled from real, binned short-read
contigs (fragmented_database_final/, one allDistances-style entry per
plasmid pair) against the ground-truth distance computed directly from the
full-length plasmid_references/ sequences for the same pair. For each pair,
takes the minimum of the binned distribution and the minimum of the
ground-truth distribution, then reports how often/by how much the binned
minimum over- or under-estimates the reference minimum, plus a scatter plot
(ground truth vs binned min) and per-pair violin plots for the mismatched
pairs.

CLI: python plot_vs_ground_trhuth.py <binned_file> <reference_file>
"""
import sys
import re
import matplotlib.pyplot as plt
import numpy as np

# ── configurable prefixes ────────────────────────────────────────────────────
FILE1_PREFIX = "/fragmented_database_final/"
FILE2_PREFIX = "/plasmid_references/"


def parse_file1(filepath):
    """
    File 1: both paths share FILE1_PREFIX, first path always ends with _binned.fasta
    Line 1: /fragmented_database_final/A_binned.fasta/fragmented_database_final/B_binned.fasta
    Line 2: space-separated floats
    Returns: { (name_A, name_B): [float, ...] }
    """
    records = {}
    with open(filepath) as f:
        lines = [l.rstrip("\n") for l in f if l.strip()]

    if len(lines) % 2 != 0:
        raise ValueError(f"{filepath}: expected even number of non-empty lines, got {len(lines)}")

    for i in range(0, len(lines), 2):
        pair_line = lines[i]
        dist_line = lines[i + 1]

        # Split on the boundary: "_binned.fasta" + FILE1_PREFIX
        # e.g. ".../A_binned.fasta/fragmented_database_final/B_binned.fasta"
        #                         ^--- split here
        split_token = f"_binned.fasta{FILE1_PREFIX}"
        idx = pair_line.find(split_token)
        if idx == -1:
            print(f"[WARN] Cannot split pair line at record {i//2 + 1}: {pair_line!r} — skipping")
            continue
        
        name1 = pair_line.split("/")[2].replace("_binned.fasta","")
        name2 = pair_line.split("/")[4].replace("_binned.fasta","")
        

        try:
            values = [float(x) for x in dist_line.split()]
        except ValueError:
            print(f"[WARN] Cannot parse distribution at record {i//2 + 1}: {dist_line!r} — skipping")
            continue

        records[(name1, name2)] = values

    return records


def parse_file2(filepath):
    """
    File 2: both paths share FILE2_PREFIX, paths end with .fasta (no _binned suffix)
    Line 1: /plasmid_references/A.fasta/plasmid_references/B.fasta
    Line 2: space-separated floats
    Returns: { (name_A, name_B): [float, ...] }
    """
    records = {}
    with open(filepath) as f:
        lines = [l.rstrip("\n") for l in f if l.strip()]

    if len(lines) % 2 != 0:
        raise ValueError(f"{filepath}: expected even number of non-empty lines, got {len(lines)}")

    for i in range(0, len(lines), 2):
        pair_line = lines[i]
        dist_line = lines[i + 1]

        # Split on ".fasta" + FILE2_PREFIX boundary
        split_token = f".fasta{FILE2_PREFIX}"
        idx = pair_line.find(split_token)
        if idx == -1:
            print(f"[WARN] Cannot split pair line at record {i//2 + 1}: {pair_line!r} — skipping")
            continue

        name1 = pair_line.split("/")[2].replace(".fasta","")
        name2 = pair_line.split("/")[4].replace(".fasta","")
        
        try:
            values = [float(x) for x in dist_line.split()]
        except ValueError:
            print(f"[WARN] Cannot parse distribution at record {i//2 + 1}: {dist_line!r} — skipping")
            continue

        records[(name1, name2)] = values

    return records


def match_records(file1_records, file2_records):
    """
    Match file1 pairs to file2 pairs by stripping '_binned' from file1 names.
    Returns: [(label, min_binned, ground_truth_min)]
    """
    matched = []

    for (n1, n2), dist in file1_records.items():
        # Normalise: CAV1668_all_references_0_binned.fasta -> CAV1668_all_references_0.fasta

        ref_key = (n1, n2)

        if ref_key not in file2_records:
            print(f"[WARN] No ground-truth match for ({n1!r}, {n2!r}) — skipping")
            continue

        min_binned   = min(dist)
        ground_truth = min(file2_records[ref_key])

        # Short readable label: "CAV1668_0 vs KPNIH27_0"
        label = f"{n1.replace('_all_references', '')} vs " \
                f"{n2.replace('_all_references', '')}"

        matched.append((label, min_binned, ground_truth,dist))

    return matched


def plot_line(matched, out_path="comparison.png"):
    """Line-plot variant (ground truth vs binned min, sorted by ground
    truth) of the comparison. Not called by __main__ (which uses plot());
    kept as an alternate visualisation."""
    if not matched:
        print("[ERROR] Nothing to plot.")
        return

    # Sort by ground truth value for visual clarity
    matched.sort(key=lambda x: x[2])

    labels    = [m[0] for m in matched]
    min_vals  = np.array([m[1] for m in matched])
    gt_vals   = np.array([m[2] for m in matched])
    x         = np.arange(len(labels))

    fig, ax = plt.subplots(figsize=(max(10, len(labels) * 0.7), 5))

    ax.plot(x, gt_vals,  "o-",  color="steelblue", linewidth=1.5, markersize=5, label="Ground truth (min)")
    ax.plot(x, min_vals, "s--", color="tomato",    linewidth=1.5, markersize=5, label="Binned min")

    ax.set_xticks(x)
    # ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("Distance")
    ax.set_xlabel("Plasmid pair")
    ax.set_title("Minimum binned distance vs ground-truth distance per plasmid pair")
    ax.legend()
    ax.grid(axis="y", linestyle="--", alpha=0.4)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    print(f"[INFO] Plot saved to {out_path}")
    plt.show()

def plot(matched, out_path_scatter="comparison_scatter.png", out_path_violin="comparison_violin_filter.pdf"):
    """Main comparison output, called from __main__.
    1. Prints summary stats: how many pairs the binned min exactly matches,
       under-, or over-estimates the ground-truth min, plus average/max
       absolute error.
    2. Scatter plot (ground truth vs binned min, y=x line for a perfect
       match; point size/annotation reflects how many pairs share the same
       (x, y) coordinate).
    3. Multi-page PDF of violin plots (one per mismatched pair), showing
       the full binned distance distribution with the binned min and
       ground truth marked, so individual mismatches can be inspected."""
    if not matched:
        print("[ERROR] Nothing to plot.")
        return

    labels   = [m[0] for m in matched]
    min_vals = np.array([m[1] for m in matched])
    gt_vals  = np.array([m[2] for m in matched])
    dists    = [m[3] for m in matched]

    # ── STATISTICS ───────────────────────────────────────────────────────────
    n_total   = len(matched)
    n_equal   = sum(mv == gv for mv, gv in zip(min_vals, gt_vals))
    n_lower   = sum(mv <  gv for mv, gv in zip(min_vals, gt_vals))
    n_higher  = sum(mv >  gv for mv, gv in zip(min_vals, gt_vals))

    diffs     = np.abs(min_vals - gt_vals)
    avg_dist  = np.mean(diffs)
    max_dist  = np.max(diffs)
    max_label = labels[np.argmax(diffs)]

    print("=" * 50)
    print(f"  Total pairs:                     {n_total}")
    print(f"  Equal to ground truth:           {n_equal}  ({100*n_equal/n_total:.1f}%)")
    print(f"  Min LOWER than ground truth:     {n_lower}  ({100*n_lower/n_total:.1f}%)")
    print(f"  Min HIGHER than ground truth:    {n_higher}  ({100*n_higher/n_total:.1f}%)")
    print(f"  Average distance to ground truth: {avg_dist:.4f}")
    print(f"  Maximum distance to ground truth: {max_dist:.4f}  ({max_label})")
    print("=" * 50)

    # ── 1. SCATTER PLOT ──────────────────────────────────────────────────────
    from collections import Counter
    coord_counts = Counter(zip(gt_vals.tolist(), min_vals.tolist()))
    sizes        = np.array([coord_counts[(gx, mx)] for gx, mx in zip(gt_vals.tolist(), min_vals.tolist())])
    point_sizes  = 40 + sizes * 60

    fig1, ax1 = plt.subplots(figsize=(7, 7))

    ax1.scatter(gt_vals, min_vals, s=point_sizes, alpha=0.6,
                edgecolors="k", linewidths=0.5, color="steelblue", zorder=3)

    for (gx, mx), count in coord_counts.items():
        if count > 1:
            ax1.annotate(f"n={count}", (gx, mx),
                         textcoords="offset points", xytext=(6, 4),
                         fontsize=7, color="navy")

    all_vals = np.concatenate([gt_vals, min_vals])
    lim_min  = all_vals.min() - 0.05 * np.ptp(all_vals)
    lim_max  = all_vals.max() + 0.05 * np.ptp(all_vals)
    ax1.plot([lim_min, lim_max], [lim_min, lim_max],
             color="tomato", linewidth=1.2, linestyle="--", label="y = x (perfect match)")

    # for label, gx, mx, sz in zip(labels, gt_vals, min_vals, sizes):
    #     if sz == 1:
    #         ax1.annotate(label, (gx, mx), textcoords="offset points",
    #                      xytext=(5, 3), fontsize=6, alpha=0.75)

    # Add stats box to scatter plot
    stats_text = (f"n={n_total}  |  equal: {n_equal} ({100*n_equal/n_total:.1f}%)  |  "
                  f"lower: {n_lower} ({100*n_lower/n_total:.1f}%)  |  "
                  f"higher: {n_higher} ({100*n_higher/n_total:.1f}%)\n"
                  f"avg Δ={avg_dist:.4f}  |  max Δ={max_dist:.4f} ({max_label})")
    ax1.text(0.02, 0.98, stats_text, transform=ax1.transAxes,
             fontsize=7, verticalalignment="top",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))

    ax1.set_xlim(lim_min, lim_max)
    ax1.set_ylim(lim_min, lim_max)
    ax1.set_xlabel("Ground truth (min distance)")
    ax1.set_ylabel("Binned min distance")
    ax1.set_title("Binned min vs ground-truth distance\n(point size reflects overlap count)")
    ax1.legend()
    ax1.grid(linestyle="--", alpha=0.3)
    ax1.set_aspect("equal")

    plt.tight_layout()
    plt.savefig(out_path_scatter, dpi=150)
    print(f"[INFO] Scatter plot saved to {out_path_scatter}")
    plt.show()
    plt.close(fig1)

    # ── 2. VIOLIN PLOTS (PDF) for unequal pairs ──────────────────────────────
    unequal = [(lab, dist, mv, gv)
               for lab, mv, gv, dist in zip(labels, min_vals, gt_vals, dists)
               if mv != gv]

    if not unequal:
        print("[INFO] All pairs equal ground truth — no violin plot needed.")
        return

    from matplotlib.backends.backend_pdf import PdfPages

    n_cols       = 3
    n_rows_page  = 3                          # subplots per PDF page
    plots_per_page = n_cols * n_rows_page

    with PdfPages(out_path_violin) as pdf:

        # ── summary page ─────────────────────────────────────────────────────
        fig_summary, ax_s = plt.subplots(figsize=(8, 3))
        ax_s.axis("off")
        summary_lines = [
            f"Total pairs:                  {n_total}",
            f"Equal to ground truth:        {n_equal}  ({100*n_equal/n_total:.1f}%)",
            f"Min LOWER than ground truth:  {n_lower}  ({100*n_lower/n_total:.1f}%)",
            f"Min HIGHER than ground truth: {n_higher}  ({100*n_higher/n_total:.1f}%)",
            f"Average distance to GT:       {avg_dist:.4f}",
            f"Maximum distance to GT:       {max_dist:.4f}  ({max_label})",
            f"",
            f"Pages below: violin distributions for the {len(unequal)} unequal pairs.",
            f"  red dot      = binned min",
            f"  green diamond = ground truth",
        ]
        ax_s.text(0.05, 0.95, "\n".join(summary_lines),
                  transform=ax_s.transAxes, fontsize=11,
                  verticalalignment="top", fontfamily="monospace",
                  bbox=dict(boxstyle="round,pad=0.5", facecolor="lightyellow", alpha=0.8))
        ax_s.set_title("Summary statistics", fontsize=13, fontweight="bold")
        pdf.savefig(fig_summary, bbox_inches="tight")
        plt.close(fig_summary)

        # ── violin pages ─────────────────────────────────────────────────────
        for page_start in range(0, len(unequal), plots_per_page):
            page_items = unequal[page_start: page_start + plots_per_page]
            n_rows     = int(np.ceil(len(page_items) / n_cols))

            fig, axes = plt.subplots(n_rows, n_cols,
                                     figsize=(n_cols * 4, n_rows * 4),
                                     squeeze=False)

            for idx, (lab, dist, min_val, gt_val) in enumerate(page_items):
                row, col = divmod(idx, n_cols)
                ax = axes[row][col]

                vp = ax.violinplot(dist, positions=[0], showmedians=True,
                                   showextrema=True, widths=0.6)

                for pc in vp["bodies"]:
                    pc.set_facecolor("steelblue")
                    pc.set_alpha(0.6)
                vp["cmedians"].set_color("navy")
                vp["cbars"].set_color("steelblue")
                vp["cmaxes"].set_color("steelblue")
                vp["cmins"].set_color("steelblue")

                ax.scatter([0], [min_val], color="tomato",    zorder=5, s=60,
                           label=f"Binned min ({min_val:.3f})")
                ax.scatter([0], [gt_val],  color="limegreen", zorder=5, s=60,
                           marker="D", label=f"Ground truth ({gt_val:.3f})")

                # Annotate whether binned min is lower or higher
                direction = "↓ lower" if min_val < gt_val else "↑ higher"
                diff      = abs(min_val - gt_val)
                ax.set_title(f"{lab}\n{direction}  |  Δ={diff:.3f}", fontsize=8)
                ax.set_xticks([])
                ax.set_ylabel("Distance", fontsize=8)
                ax.legend(fontsize=6)
                ax.grid(axis="y", linestyle="--", alpha=0.3)

            # Hide unused subplots on last page
            for idx in range(len(page_items), n_rows * n_cols):
                row, col = divmod(idx, n_cols)
                axes[row][col].set_visible(False)

            fig.suptitle(
                f"Unequal pairs — page {page_start // plots_per_page + 1} of "
                f"{int(np.ceil(len(unequal) / plots_per_page))}",
                fontsize=11)
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    print(f"[INFO] Violin PDF saved to {out_path_violin}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: python {sys.argv[0]} <binned_file> <reference_file>")
        sys.exit(1)

    print("[INFO] Parsing binned file …")
    file1_records = parse_file1(sys.argv[1])
    print(f"       {len(file1_records)} records loaded")

    print("[INFO] Parsing reference file …")
    file2_records = parse_file2(sys.argv[2])
    print(f"       {len(file2_records)} records loaded")

    matched = match_records(file1_records, file2_records)
    print(f"[INFO] {len(matched)} matched pairs — plotting …")

    plot(matched)