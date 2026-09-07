
"""
Parses a binned fragments file where pairs of FASTA paths are followed by a
shared line of bin numbers. Each pair is treated as one plasmid pair data point.
Computes the total sequence length of each FASTA file from two different folders,
averages across the two files in each pair per folder, then takes the absolute
difference between the two folders. Produces a scatter plot of that difference (x)
vs bin number variance (y), one point per plasmid pair, with Spearman rank
correlation and p-value.

Also reports, as a table and as an average, how much sequence is "missing" per
file — i.e. the absolute difference in total length between a file's version
in folder A and its version in folder B.

Expected input file structure (two paths, then their shared numbers, repeating):
    /binned_fragments/CAV1321_all_references_6_binned.fasta/binned_fragments/CAV1344_all_references_3_binned.fasta
    6 6 4 7 7 6 6 6 ...
"""

import sys
import os
import re
import csv
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


def find_fasta_in_db(fasta_path, db_dir, kind):
    """Walk db_dir to find a file matching the basename of fasta_path."""
    basename = os.path.basename(fasta_path)
    if kind == "a":
        basename = basename.replace("_binned", "")
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


def avg_total_length_for_pair(fasta_paths, db_dir, kind):
    """
    Compute average total sequence length across all files in a pair,
    looked up in db_dir. Returns None if any file is missing.
    """
    lengths = []
    for fp in fasta_paths:
        full_path = find_fasta_in_db(fp, db_dir, kind)
        if full_path is None:
            return None
        lengths.append(total_sequence_length(full_path))
    return float(np.mean(lengths))


def pair_label(paths):
    """Build a short readable label from two FASTA basenames."""
    ids = []
    for p in paths:
        m = re.search(r'(CAV\d+)', os.path.basename(p))
        if m:
            ids.append(m.group(1))
    return " / ".join(ids) if ids else os.path.basename(paths[0])


# ── per-file length-difference ("missing") report ───────────────────────────

def per_file_length_report(entries, folder_a, folder_b):
    """
    For every individual FASTA file referenced in the binned file (deduplicated),
    look up its total sequence length in folder_a and in folder_b, then compute
    the absolute difference between the two ("how much is missing" for that file).

    Returns a list of dicts with keys: file, length_a, length_b, missing_bp,
    pct_missing, note. length_a/length_b/missing_bp/pct_missing are None if the
    file couldn't be found in one or both folders.
    """
    seen = set()
    rows = []

    for fasta_paths, _ in entries:
        for fp in fasta_paths:
            if fp in seen:
                continue
            seen.add(fp)

            path_a = find_fasta_in_db(fp, folder_a, "a")
            path_b = find_fasta_in_db(fp, folder_b, "b")
            basename = os.path.basename(fp)

            if path_a is None or path_b is None:
                if path_a is None and path_b is None:
                    note = "not found in either folder"
                elif path_a is None:
                    note = "not found in folder A"
                else:
                    note = "not found in folder B"
                rows.append({
                    "file": basename, "length_a": None, "length_b": None,
                    "missing_bp": None, "pct_missing": None, "note": note,
                })
                continue

            length_a = total_sequence_length(path_a)
            length_b = total_sequence_length(path_b)
            missing_bp = abs(length_a - length_b)
            pct_missing = (missing_bp / length_a * 100) if length_a else 0.0

            rows.append({
                "file": basename, "length_a": length_a, "length_b": length_b,
                "missing_bp": missing_bp, "pct_missing": pct_missing, "note": "",
            })

    return rows


def print_length_diff_table(rows, folder_a, folder_b):
    """
    Prints a per-file table of length differences ("missing" bp) between
    folder_a and folder_b, followed by the average across all files that had
    data in both folders. Returns the summary stats as a dict.
    """
    if not rows:
        print("No files to report on.")
        return None

    a_name = os.path.basename(folder_a.rstrip("/")) or folder_a
    b_name = os.path.basename(folder_b.rstrip("/")) or folder_b

    headers = ["File", f"Length in {a_name} (bp)", f"Length in {b_name} (bp)", "Missing (bp)", "% Missing"]
    col_widths = [
        max(len(headers[0]), max(len(r["file"]) for r in rows)),
        len(headers[1]),
        len(headers[2]),
        len(headers[3]),
        len(headers[4]),
    ]

    def fmt_row(values):
        return "  ".join(str(v).ljust(w) for v, w in zip(values, col_widths))

    width = sum(col_widths) + 2 * (len(col_widths) - 1)

    print("\n" + "=" * width)
    print("Per-file length difference (\"missing\" sequence) report")
    print("=" * width)
    print(fmt_row(headers))
    print("-" * width)

    for r in rows:
        if r["length_a"] is None:
            print(fmt_row([r["file"], "N/A", "N/A", "N/A", r["note"]]))
        else:
            print(fmt_row([
                r["file"],
                f"{r['length_a']:,}",
                f"{r['length_b']:,}",
                f"{r['missing_bp']:,}",
                f"{r['pct_missing']:.2f}%",
            ]))

    valid = [r for r in rows if r["missing_bp"] is not None]
    print("-" * width)

    if not valid:
        print("No files had length data in both folders — cannot compute an average.")
        print("=" * width + "\n")
        return None

    avg_missing_bp = float(np.mean([r["missing_bp"] for r in valid]))
    avg_pct_missing = float(np.mean([r["pct_missing"] for r in valid]))
    total_missing_bp = sum(r["missing_bp"] for r in valid)
    total_length_a = sum(r["length_a"] for r in valid)
    overall_pct_missing = (total_missing_bp / total_length_a * 100) if total_length_a else 0.0

    print(fmt_row(["AVERAGE", "", "", f"{avg_missing_bp:,.1f}", f"{avg_pct_missing:.2f}%"]))
    print("=" * width)
    print(f"Across {len(valid)} file(s) with data in both folders: "
          f"average missing = {avg_missing_bp:,.1f} bp ({avg_pct_missing:.2f}% of length in {a_name} on average); "
          f"aggregate missing = {total_missing_bp:,} bp ({overall_pct_missing:.2f}% of total length in {a_name}).\n")

    if len(valid) < len(rows):
        print(f"Note: {len(rows) - len(valid)} file(s) were excluded from the average "
              f"because they were not found in both folders.\n")

    return {
        "avg_missing_bp": avg_missing_bp,
        "avg_pct_missing": avg_pct_missing,
        "total_missing_bp": total_missing_bp,
        "overall_pct_missing": overall_pct_missing,
        "n_files": len(valid),
    }


def save_length_diff_csv(rows, out_path):
    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["file", "length_a_bp", "length_b_bp", "missing_bp", "pct_missing", "note"])
        for r in rows:
            writer.writerow([
                r["file"],
                r["length_a"] if r["length_a"] is not None else "",
                r["length_b"] if r["length_b"] is not None else "",
                r["missing_bp"] if r["missing_bp"] is not None else "",
                f"{r['pct_missing']:.4f}" if r["pct_missing"] is not None else "",
                r["note"],
            ])
    print(f"Per-file length-difference report saved to: {out_path}")


# ── plotting ──────────────────────────────────────────────────────────────────

def plot_scatter(points, folder_a, folder_b):
    """
    points: list of dicts with keys 'label', 'length_diff', 'variance'
    """
    x = np.array([p["length_diff"] for p in points])
    y = np.array([p["variance"] for p in points])
    labels = [p["label"] for p in points]

    fig, ax = plt.subplots(figsize=(9, 6))

    scatter = ax.scatter(
        x, y,
        s=90, zorder=3,
        c=y, cmap="viridis",
        edgecolors="white", linewidths=0.6
    )

    cbar = fig.colorbar(scatter, ax=ax, pad=0.02)
    cbar.set_label("Variance", fontsize=9)
    cbar.ax.yaxis.set_tick_params()
    plt.setp(cbar.ax.yaxis.get_ticklabels())

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
        ax.text(
            0.05, 0.95, annotation,
            transform=ax.transAxes,
            fontsize=9, color="#ffa657",
            verticalalignment="top",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#1a1d27", edgecolor="#ffa657", alpha=0.8)
        )
        print(f"Spearman rank correlation (ρ): {r:.4f}")
        print(f"p-value: {p_value:.4f} ({sig_label} at α=0.05)")
        print(f"Note: based on subsample of {len(x)} pairs — interpret p-value with caution.")
    elif len(x) == 2:
        r, _ = stats.spearmanr(x, y)
        ax.text(
            0.05, 0.95, f"Spearman ρ = {float(r):.4f}\n(need ≥3 pairs for p-value)",
            transform=ax.transAxes,
            fontsize=9, color="#ffa657",
            verticalalignment="top",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#1a1d27", edgecolor="#ffa657", alpha=0.8)
        )
        print(f"Spearman ρ = {float(r):.4f} (need at least 3 pairs to compute p-value)")
    else:
        print("Not enough points to compute correlation (need at least 2).")

    a_name = os.path.basename(folder_a.rstrip("/"))
    b_name = os.path.basename(folder_b.rstrip("/"))
    ax.set_xlabel(f"Absolute difference in avg total length (bp)\n|{a_name} − {b_name}|", fontsize=11)
    ax.set_ylabel("Variance of bin number distribution", fontsize=11)
    ax.set_title("Plasmid pair: length difference between folders vs bin variance", fontsize=13)
    ax.grid(True, color="#30363d", linewidth=0.5, linestyle="--")

    plt.tight_layout()
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "binned_analysis.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    print(f"Plot saved to: {out_path}")
    plt.show()


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) != 4:
        print("Usage: python analyze_binned_fragments.py <binned_file> <folder_a> <folder_b>")
        sys.exit(1)

    binned_file = sys.argv[1]
    folder_a = sys.argv[2]
    folder_b = sys.argv[3]

    for path, label in [(binned_file, "binned file"), (folder_a, "folder A"), (folder_b, "folder B")]:
        if not os.path.exists(path):
            print(f"Error: {label} not found: {path}")
            sys.exit(1)

    print(f"Parsing: {binned_file}")
    entries = parse_binned_file(binned_file)
    print(f"Found {len(entries)} plasmid pair(s).\n")

    # Per-file "missing" (length difference) report: table + averages
    file_rows = per_file_length_report(entries, folder_a, folder_b)
    print_length_diff_table(file_rows, folder_a, folder_b)
    length_csv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "length_diff_report.csv")
    save_length_diff_csv(file_rows, length_csv_path)

    points = []
    for fasta_paths, numbers in entries:
        label = pair_label(fasta_paths)

        avg_a = avg_total_length_for_pair(fasta_paths, folder_a, "a")
        if avg_a is None:
            print(f"  WARNING: pair [{label}] not fully found in folder A. Skipping from plot.")
            continue

        avg_b = avg_total_length_for_pair(fasta_paths, folder_b, "b")
        if avg_b is None:
            print(f"  WARNING: pair [{label}] not fully found in folder B. Skipping from plot.")
            continue

        length_diff = abs(avg_a - avg_b)
        variance = float(np.var(numbers))
        print(f"  Pair [{label}]: avg length A = {avg_a:,.1f} bp, avg length B = {avg_b:,.1f} bp, "
              f"|diff| = {length_diff:,.1f} bp, variance = {variance:.4f}")
        points.append({"label": label, "length_diff": length_diff, "variance": variance})

    if not points:
        print("No valid plasmid pairs found (all had missing files). Cannot plot.")
        sys.exit(1)

    print(f"\nPlotting {len(points)} pairs...")
    plot_scatter(points, folder_a, folder_b)


if __name__ == "__main__":
    main()