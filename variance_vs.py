"""
variance_vs.py — variance of the binned distance distribution vs various
plasmid-pair properties.

This file is a stack of six near-identical script variants, each exploring
whether the variance of a pair's binned short-read distance distribution
(the spread in allDistances_v2_*.txt-style output, i.e. how sensitive the
result is to fragment ordering) correlates with some property of the pair:
number of contigs, total plasmid length, length difference between two
database versions, mash distance, and (further down) the *error* against
ground truth rather than variance itself, plus a standalone pling-vs-mash
scatter plot. Each variant is a full commented-out module (its own
docstring, argument parsing, plotting) stacked one after another; only the
"variance vs total length" block (the one active section, starting at
`## This is for variance vs total length` below) is uncommented and
actually runs via `if __name__ == "__main__"`. To use one of the other
variants, uncomment that block and comment out (or otherwise disable) the
active one, since they'd otherwise both try to define `main()`.
"""

## This is for variance vs number of fragments

# """
# analyze_binned_fragments.py
 
# Parses a binned fragments file where pairs of FASTA paths are followed by a
# shared line of bin numbers. Each pair is treated as one plasmid pair data point.
# Counts '>' header lines in each FASTA file found under 'fragmented_database_final',
# averages the two counts, then produces a scatter plot of average header count (x)
# vs bin number variance (y), one point per plasmid pair.
 
# Expected input file structure (two paths, then their shared numbers, repeating):
#     /binned_fragments/CAV1321_all_references_6_binned.fasta/binned_fragments/CAV1344_all_references_3_binned.fasta
#     6 6 4 7 7 6 6 6 ...
#     /binned_fragments/CAV1400_all_references_2_binned.fasta/binned_fragments/CAV1401_all_references_5_binned.fasta
#     5 5 6 7 6 5 7 6 ...
 
# Usage:
#     python analyze_binned_fragments.py <binned_file> <fragmented_database_final_dir>
 
# Example:
#     python analyze_binned_fragments.py my_binned.txt fragmented_database_final/
# """
 
# import sys
# import os
# import re
# import numpy as np
# import matplotlib.pyplot as plt
# from scipy import stats
 
 
# # ── helpers ───────────────────────────────────────────────────────────────────
 
# def parse_binned_file(filepath):
#     """
#     Parse the binned file into a list of ([fasta_paths], [numbers]) tuples.
#     Each group of FASTA paths (possibly concatenated on one line) is paired
#     with the number sequence on the following line.
#     """
#     entries = []
#     current_paths = []
 
#     with open(filepath) as fh:
#         for line in fh:
#             line = line.strip()
#             if not line:
#                 continue
 
#             if ".fasta" in line:
#                 # Split concatenated paths like /a/b.fasta/c/d.fasta
#                 parts = re.split(r'(?=/)', line)
#                 for part in parts:
#                     part = part.strip()
#                     if part.endswith(".fasta"):
#                         current_paths.append(part)
#             else:
#                 tokens = line.split()
#                 numbers = []
#                 for t in tokens:
#                     try:
#                         numbers.append(int(t))
#                     except ValueError:
#                         pass
 
#                 if numbers and current_paths:
#                     entries.append((list(current_paths), numbers))
#                     current_paths = []
 
#     return entries
 
 
# def find_fasta_in_db(fasta_path, db_dir):
#     """Walk db_dir to find a file matching the basename of fasta_path."""
#     basename = os.path.basename(fasta_path)
#     for root, dirs, files in os.walk(db_dir):
#         if basename in files:
#             return os.path.join(root, basename)
#     return None
 
 
# def count_fasta_headers(filepath):
#     """Count lines starting with '>' in a FASTA file."""
#     count = 0
#     with open(filepath) as fh:
#         for line in fh:
#             if line.startswith(">"):
#                 count += 1
#     return count
 
 
# def pair_label(paths):
#     """Build a short readable label from two FASTA basenames."""
#     ids = []
#     for p in paths:
#         m = re.search(r'(CAV\d+)', os.path.basename(p))
#         if m:
#             ids.append(m.group(1))
#     return " / ".join(ids) if ids else os.path.basename(paths[0])
 
 
# # ── plotting ──────────────────────────────────────────────────────────────────
 
# def plot_scatter(points):
#     """
#     points: list of dicts with keys 'label', 'avg_header_count', 'variance'
#     """
#     x = np.array([p["avg_header_count"] for p in points])
#     y = np.array([p["variance"] for p in points])
#     labels = [p["label"] for p in points]

#     fig, ax = plt.subplots(figsize=(9, 6))

#     scatter = ax.scatter(
#         x, y,
#         s=90, zorder=3,
#         edgecolors="black", linewidths=0.6
#     )

#     # cbar = fig.colorbar(scatter, ax=ax, pad=0.02)
#     # cbar.set_label("Variance", fontsize=9)

#     # Spearman rank correlation + significance test
#     if len(x) >= 3:
#         r, p_value = stats.spearmanr(x, y)
#         r, p_value = float(r), float(p_value)
#         sig_label = "significant" if p_value < 0.05 else "not significant"
#         annotation = (
#             f"Spearman ρ = {r:.4f}\n"
#             f"p-value = {p_value:.4f}\n"
#             f"({sig_label} at α=0.05)\n"
#             f"n = {len(x)} pairs (subsample)"
#         )
#         # ax.text(
#         #     0.05, 0.95, annotation,
#         #     transform=ax.transAxes,
#         #     fontsize=9, color="black",
#         #     verticalalignment="top",
#         #     bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="black", alpha=0.8)
#         # )
#         print(f"Spearman rank correlation (ρ): {r:.4f}")
#         print(f"p-value: {p_value:.4f} ({sig_label} at α=0.05)")
#         print(f"Note: based on subsample of {len(x)} pairs — interpret p-value with caution.")
#     elif len(x) == 2:
#         r, _ = stats.spearmanr(x, y)
#         # ax.text(
#         #     0.05, 0.95, f"Spearman ρ = {float(r):.4f}\n(need ≥3 pairs for p-value)",
#         #     transform=ax.transAxes,
#         #     fontsize=9, color="black",
#         #     verticalalignment="top",
#         #     bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="black", alpha=0.8)
#         # )
#         print(f"Spearman ρ = {float(r):.4f} (need at least 3 pairs to compute p-value)")
#     else:
#         print("Not enough points to compute correlation (need at least 2).")

#     ax.set_xlabel("Average number of contigs", fontsize=11)
#     ax.set_ylabel("Variance of Pling distances", fontsize=11)
#     ax.set_title("", fontsize=13)
#     ax.grid(True, color="lightgrey", linewidth=0.5, linestyle="--")

#     plt.tight_layout()
#     out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "binned_analysis.png")
#     plt.savefig(out_path, dpi=150, bbox_inches="tight")
#     print(f"Plot saved to: {out_path}")
#     plt.show()
 
 
# # ── main ──────────────────────────────────────────────────────────────────────
 
# def main():
#     if len(sys.argv) != 3:
#         print("Usage: python analyze_binned_fragments.py <binned_file> <fragmented_database_final_dir>")
#         sys.exit(1)
 
#     binned_file = sys.argv[1]
#     db_dir = sys.argv[2]
 
#     if not os.path.isfile(binned_file):
#         print(f"Error: binned file not found: {binned_file}")
#         sys.exit(1)
#     if not os.path.isdir(db_dir):
#         print(f"Error: database directory not found: {db_dir}")
#         sys.exit(1)
 
#     print(f"Parsing: {binned_file}")
#     entries = parse_binned_file(binned_file)
#     print(f"Found {len(entries)} plasmid pair(s).")
 
#     points = []
#     for fasta_paths, numbers in entries:
#         label = pair_label(fasta_paths)
#         header_counts = []
 
#         for fp in fasta_paths:
#             name = os.path.basename(fp)
#             full_path = find_fasta_in_db(fp, db_dir)
#             if full_path is None:
#                 print(f"  WARNING: '{name}' not found in '{db_dir}'. Skipping pair.")
#                 break
#             hcount = count_fasta_headers(full_path)
#             print(f"  {name}: {hcount} sequences")
#             header_counts.append(hcount)
#         else:
#             avg_count = float(np.mean(header_counts))
#             variance = float(np.var(numbers))
#             print(f"  Pair [{label}]: avg headers = {avg_count:.1f}, variance = {variance:.4f} (n={len(numbers)} bins)\n")
#             points.append({"label": label, "avg_header_count": avg_count, "variance": variance})
 
#     if not points:
#         print("No valid plasmid pairs found. Cannot plot.")
#         sys.exit(1)
 
#     plot_scatter(points)
 
 
# if __name__ == "__main__":
#     main()

## This is for variance vs total length

"""
analyze_binned_fragments.py
 
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
 
Usage:
    python analyze_binned_fragments.py <binned_file> <fragmented_database_final_dir>
 
Example:
    python analyze_binned_fragments.py my_binned.txt fragmented_database_final/
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
 
## This is for variance vs the difference between the two

# """
# analyze_binned_fragments.py

# Parses a binned fragments file where pairs of FASTA paths are followed by a
# shared line of bin numbers. Each pair is treated as one plasmid pair data point.
# Computes the total sequence length of each FASTA file from two different folders,
# averages across the two files in each pair per folder, then takes the absolute
# difference between the two folders. Produces a scatter plot of that difference (x)
# vs bin number variance (y), one point per plasmid pair, with Spearman rank
# correlation and p-value.

# Also reports, as a table and as an average, how much sequence is "missing" per
# file — i.e. the absolute difference in total length between a file's version
# in folder A and its version in folder B.

# Expected input file structure (two paths, then their shared numbers, repeating):
#     /binned_fragments/CAV1321_all_references_6_binned.fasta/binned_fragments/CAV1344_all_references_3_binned.fasta
#     6 6 4 7 7 6 6 6 ...

# Usage:
#     python analyze_binned_fragments.py <binned_file> <folder_a> <folder_b>

# Example:
#     python analyze_binned_fragments.py my_binned.txt fragmented_database_final/ fragmented_database_alt/
# """

# import sys
# import os
# import re
# import csv
# import numpy as np
# import matplotlib.pyplot as plt
# from scipy import stats


# # ── helpers ───────────────────────────────────────────────────────────────────

# def parse_binned_file(filepath):
#     entries = []
#     current_paths = []

#     with open(filepath) as fh:
#         for line in fh:
#             line = line.strip()
#             if not line:
#                 continue

#             if ".fasta" in line:
#                 parts = re.split(r'(?=/)', line)
#                 for part in parts:
#                     part = part.strip()
#                     if part.endswith(".fasta"):
#                         current_paths.append(part)
#             else:
#                 tokens = line.split()
#                 numbers = []
#                 for t in tokens:
#                     try:
#                         numbers.append(int(t))
#                     except ValueError:
#                         pass

#                 if numbers and current_paths:
#                     entries.append((list(current_paths), numbers))
#                     current_paths = []

#     return entries


# def find_fasta_in_db(fasta_path, db_dir, kind):
#     """Walk db_dir to find a file matching the basename of fasta_path."""
#     basename = os.path.basename(fasta_path)
#     if kind == "a":
#         basename = basename.replace("_binned", "")
#     for root, dirs, files in os.walk(db_dir):
#         if basename in files:
#             return os.path.join(root, basename)
#     return None


# def total_sequence_length(filepath):
#     """Sum the lengths of all sequences in a multi-FASTA file (ignores header lines)."""
#     total = 0
#     with open(filepath) as fh:
#         for line in fh:
#             line = line.strip()
#             if line and not line.startswith(">"):
#                 total += len(line)
#     return total


# def avg_total_length_for_pair(fasta_paths, db_dir, kind):
#     """
#     Compute average total sequence length across all files in a pair,
#     looked up in db_dir. Returns None if any file is missing.
#     """
#     lengths = []
#     for fp in fasta_paths:
#         full_path = find_fasta_in_db(fp, db_dir, kind)
#         if full_path is None:
#             return None
#         lengths.append(total_sequence_length(full_path))
#     return float(np.mean(lengths))


# def pair_label(paths):
#     """Build a short readable label from two FASTA basenames."""
#     ids = []
#     for p in paths:
#         m = re.search(r'(CAV\d+)', os.path.basename(p))
#         if m:
#             ids.append(m.group(1))
#     return " / ".join(ids) if ids else os.path.basename(paths[0])


# # ── per-file length-difference ("missing") report ───────────────────────────

# def per_file_length_report(entries, folder_a, folder_b):
#     """
#     For every individual FASTA file referenced in the binned file (deduplicated),
#     look up its total sequence length in folder_a and in folder_b, then compute
#     the absolute difference between the two ("how much is missing" for that file).

#     Returns a list of dicts with keys: file, length_a, length_b, missing_bp,
#     pct_missing, note. length_a/length_b/missing_bp/pct_missing are None if the
#     file couldn't be found in one or both folders.
#     """
#     seen = set()
#     rows = []

#     for fasta_paths, _ in entries:
#         for fp in fasta_paths:
#             if fp in seen:
#                 continue
#             seen.add(fp)

#             path_a = find_fasta_in_db(fp, folder_a, "a")
#             path_b = find_fasta_in_db(fp, folder_b, "b")
#             basename = os.path.basename(fp)

#             if path_a is None or path_b is None:
#                 if path_a is None and path_b is None:
#                     note = "not found in either folder"
#                 elif path_a is None:
#                     note = "not found in folder A"
#                 else:
#                     note = "not found in folder B"
#                 rows.append({
#                     "file": basename, "length_a": None, "length_b": None,
#                     "missing_bp": None, "pct_missing": None, "note": note,
#                 })
#                 continue

#             length_a = total_sequence_length(path_a)
#             length_b = total_sequence_length(path_b)
#             missing_bp = abs(length_a - length_b)
#             pct_missing = (missing_bp / length_a * 100) if length_a else 0.0

#             rows.append({
#                 "file": basename, "length_a": length_a, "length_b": length_b,
#                 "missing_bp": missing_bp, "pct_missing": pct_missing, "note": "",
#             })

#     return rows


# def print_length_diff_table(rows, folder_a, folder_b):
#     """
#     Prints a per-file table of length differences ("missing" bp) between
#     folder_a and folder_b, followed by the average across all files that had
#     data in both folders. Returns the summary stats as a dict.
#     """
#     if not rows:
#         print("No files to report on.")
#         return None

#     a_name = os.path.basename(folder_a.rstrip("/")) or folder_a
#     b_name = os.path.basename(folder_b.rstrip("/")) or folder_b

#     headers = ["File", f"Length in {a_name} (bp)", f"Length in {b_name} (bp)", "Missing (bp)", "% Missing"]
#     col_widths = [
#         max(len(headers[0]), max(len(r["file"]) for r in rows)),
#         len(headers[1]),
#         len(headers[2]),
#         len(headers[3]),
#         len(headers[4]),
#     ]

#     def fmt_row(values):
#         return "  ".join(str(v).ljust(w) for v, w in zip(values, col_widths))

#     width = sum(col_widths) + 2 * (len(col_widths) - 1)

#     print("\n" + "=" * width)
#     print("Per-file length difference (\"missing\" sequence) report")
#     print("=" * width)
#     print(fmt_row(headers))
#     print("-" * width)

#     for r in rows:
#         if r["length_a"] is None:
#             print(fmt_row([r["file"], "N/A", "N/A", "N/A", r["note"]]))
#         else:
#             print(fmt_row([
#                 r["file"],
#                 f"{r['length_a']:,}",
#                 f"{r['length_b']:,}",
#                 f"{r['missing_bp']:,}",
#                 f"{r['pct_missing']:.2f}%",
#             ]))

#     valid = [r for r in rows if r["missing_bp"] is not None]
#     print("-" * width)

#     if not valid:
#         print("No files had length data in both folders — cannot compute an average.")
#         print("=" * width + "\n")
#         return None

#     avg_missing_bp = float(np.mean([r["missing_bp"] for r in valid]))
#     avg_pct_missing = float(np.mean([r["pct_missing"] for r in valid]))
#     total_missing_bp = sum(r["missing_bp"] for r in valid)
#     total_length_a = sum(r["length_a"] for r in valid)
#     overall_pct_missing = (total_missing_bp / total_length_a * 100) if total_length_a else 0.0

#     print(fmt_row(["AVERAGE", "", "", f"{avg_missing_bp:,.1f}", f"{avg_pct_missing:.2f}%"]))
#     print("=" * width)
#     print(f"Across {len(valid)} file(s) with data in both folders: "
#           f"average missing = {avg_missing_bp:,.1f} bp ({avg_pct_missing:.2f}% of length in {a_name} on average); "
#           f"aggregate missing = {total_missing_bp:,} bp ({overall_pct_missing:.2f}% of total length in {a_name}).\n")

#     if len(valid) < len(rows):
#         print(f"Note: {len(rows) - len(valid)} file(s) were excluded from the average "
#               f"because they were not found in both folders.\n")

#     return {
#         "avg_missing_bp": avg_missing_bp,
#         "avg_pct_missing": avg_pct_missing,
#         "total_missing_bp": total_missing_bp,
#         "overall_pct_missing": overall_pct_missing,
#         "n_files": len(valid),
#     }


# def save_length_diff_csv(rows, out_path):
#     with open(out_path, "w", newline="") as fh:
#         writer = csv.writer(fh)
#         writer.writerow(["file", "length_a_bp", "length_b_bp", "missing_bp", "pct_missing", "note"])
#         for r in rows:
#             writer.writerow([
#                 r["file"],
#                 r["length_a"] if r["length_a"] is not None else "",
#                 r["length_b"] if r["length_b"] is not None else "",
#                 r["missing_bp"] if r["missing_bp"] is not None else "",
#                 f"{r['pct_missing']:.4f}" if r["pct_missing"] is not None else "",
#                 r["note"],
#             ])
#     print(f"Per-file length-difference report saved to: {out_path}")


# # ── plotting ──────────────────────────────────────────────────────────────────

# def plot_scatter(points, folder_a, folder_b):
#     """
#     points: list of dicts with keys 'label', 'length_diff', 'variance'
#     """
#     x = np.array([p["length_diff"] for p in points])
#     y = np.array([p["variance"] for p in points])
#     labels = [p["label"] for p in points]

#     fig, ax = plt.subplots(figsize=(9, 6))

#     scatter = ax.scatter(
#         x, y,
#         s=90, zorder=3,
#         c=y, cmap="viridis",
#         edgecolors="white", linewidths=0.6
#     )

#     cbar = fig.colorbar(scatter, ax=ax, pad=0.02)
#     cbar.set_label("Variance", fontsize=9)
#     cbar.ax.yaxis.set_tick_params()
#     plt.setp(cbar.ax.yaxis.get_ticklabels())

#     # Spearman rank correlation + significance test
#     if len(x) >= 3:
#         r, p_value = stats.spearmanr(x, y)
#         r, p_value = float(r), float(p_value)
#         sig_label = "significant" if p_value < 0.05 else "not significant"
#         annotation = (
#             f"Spearman ρ = {r:.4f}\n"
#             f"p-value = {p_value:.4f}\n"
#             f"({sig_label} at α=0.05)\n"
#             f"n = {len(x)} pairs (subsample)"
#         )
#         ax.text(
#             0.05, 0.95, annotation,
#             transform=ax.transAxes,
#             fontsize=9, color="#ffa657",
#             verticalalignment="top",
#             bbox=dict(boxstyle="round,pad=0.4", facecolor="#1a1d27", edgecolor="#ffa657", alpha=0.8)
#         )
#         print(f"Spearman rank correlation (ρ): {r:.4f}")
#         print(f"p-value: {p_value:.4f} ({sig_label} at α=0.05)")
#         print(f"Note: based on subsample of {len(x)} pairs — interpret p-value with caution.")
#     elif len(x) == 2:
#         r, _ = stats.spearmanr(x, y)
#         ax.text(
#             0.05, 0.95, f"Spearman ρ = {float(r):.4f}\n(need ≥3 pairs for p-value)",
#             transform=ax.transAxes,
#             fontsize=9, color="#ffa657",
#             verticalalignment="top",
#             bbox=dict(boxstyle="round,pad=0.4", facecolor="#1a1d27", edgecolor="#ffa657", alpha=0.8)
#         )
#         print(f"Spearman ρ = {float(r):.4f} (need at least 3 pairs to compute p-value)")
#     else:
#         print("Not enough points to compute correlation (need at least 2).")

#     a_name = os.path.basename(folder_a.rstrip("/"))
#     b_name = os.path.basename(folder_b.rstrip("/"))
#     ax.set_xlabel(f"Absolute difference in avg total length (bp)\n|{a_name} − {b_name}|", fontsize=11)
#     ax.set_ylabel("Variance of bin number distribution", fontsize=11)
#     ax.set_title("Plasmid pair: length difference between folders vs bin variance", fontsize=13)
#     ax.grid(True, color="#30363d", linewidth=0.5, linestyle="--")

#     plt.tight_layout()
#     out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "binned_analysis.png")
#     plt.savefig(out_path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
#     print(f"Plot saved to: {out_path}")
#     plt.show()


# # ── main ──────────────────────────────────────────────────────────────────────

# def main():
#     if len(sys.argv) != 4:
#         print("Usage: python analyze_binned_fragments.py <binned_file> <folder_a> <folder_b>")
#         sys.exit(1)

#     binned_file = sys.argv[1]
#     folder_a = sys.argv[2]
#     folder_b = sys.argv[3]

#     for path, label in [(binned_file, "binned file"), (folder_a, "folder A"), (folder_b, "folder B")]:
#         if not os.path.exists(path):
#             print(f"Error: {label} not found: {path}")
#             sys.exit(1)

#     print(f"Parsing: {binned_file}")
#     entries = parse_binned_file(binned_file)
#     print(f"Found {len(entries)} plasmid pair(s).\n")

#     # Per-file "missing" (length difference) report: table + averages
#     file_rows = per_file_length_report(entries, folder_a, folder_b)
#     print_length_diff_table(file_rows, folder_a, folder_b)
#     length_csv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "length_diff_report.csv")
#     save_length_diff_csv(file_rows, length_csv_path)

#     points = []
#     for fasta_paths, numbers in entries:
#         label = pair_label(fasta_paths)

#         avg_a = avg_total_length_for_pair(fasta_paths, folder_a, "a")
#         if avg_a is None:
#             print(f"  WARNING: pair [{label}] not fully found in folder A. Skipping from plot.")
#             continue

#         avg_b = avg_total_length_for_pair(fasta_paths, folder_b, "b")
#         if avg_b is None:
#             print(f"  WARNING: pair [{label}] not fully found in folder B. Skipping from plot.")
#             continue

#         length_diff = abs(avg_a - avg_b)
#         variance = float(np.var(numbers))
#         print(f"  Pair [{label}]: avg length A = {avg_a:,.1f} bp, avg length B = {avg_b:,.1f} bp, "
#               f"|diff| = {length_diff:,.1f} bp, variance = {variance:.4f}")
#         points.append({"label": label, "length_diff": length_diff, "variance": variance})

#     if not points:
#         print("No valid plasmid pairs found (all had missing files). Cannot plot.")
#         sys.exit(1)

#     print(f"\nPlotting {len(points)} pairs...")
#     plot_scatter(points, folder_a, folder_b)


# if __name__ == "__main__":
#     main()

## This is variance vs mash distance

# """
# analyze_binned_fragments.py
 
# Parses a binned fragments file where pairs of FASTA paths are followed by a
# shared line of bin numbers. Each pair is treated as one plasmid pair data point.
# Looks up the mash distance for each pair from a mash distance file (tab-separated,
# distance in column 3), matching by filename after replacing 'plasmid_references'
# with 'binned_fragments'. Produces a scatter plot of mash distance (x) vs bin
# number variance (y), one point per plasmid pair, with Spearman rank correlation
# and p-value.
 
# Expected binned file structure:
#     /binned_fragments/CAV1321_all_references_6_binned.fasta/binned_fragments/CAV1344_all_references_3_binned.fasta
#     6 6 4 7 7 6 6 6 ...
 
# Expected mash distance file structure (tab-separated):
#     ./plasmid_references/A.fasta    ./plasmid_references/B.fasta    0.012    0    950/1000
 
# Usage:
#     python analyze_binned_fragments.py <binned_file> <mash_distance_file>
 
# Example:
#     python analyze_binned_fragments.py my_binned.txt mash_distances.tsv
# """
 
# import sys
# import os
# import re
# import numpy as np
# import matplotlib.pyplot as plt
# from scipy import stats
 
 
# # ── helpers ───────────────────────────────────────────────────────────────────
 
# def parse_binned_file(filepath):
#     entries = []
#     current_paths = []
 
#     with open(filepath) as fh:
#         for line in fh:
#             line = line.strip()
#             if not line:
#                 continue
 
#             if ".fasta" in line:
#                 parts = re.split(r'(?=/)', line)
#                 for part in parts:
#                     part = part.strip()
#                     if part.endswith(".fasta"):
#                         current_paths.append(part)
#             else:
#                 tokens = line.split()
#                 numbers = []
#                 for t in tokens:
#                     try:
#                         numbers.append(int(t))
#                     except ValueError:
#                         pass
 
#                 if numbers and current_paths:
#                     entries.append((list(current_paths), numbers))
#                     current_paths = []
 
#     return entries
 
 
# def normalise_path(path):
#     """
#     Strip leading ./ and directory, return just the basename with
#     'plasmid_references' replaced by 'binned_fragments' for matching.
#     """
#     path = path.strip()
#     basename = os.path.basename(path)
#     return basename
 
 
# def parse_mash_file(filepath):
#     """
#     Parse mash distance file into a dict:
#         (basename_a, basename_b) -> mash_distance
#     Stored in both orders for easy lookup.
#     'plasmid_references' is replaced with 'binned_fragments' in basenames.
#     """
#     distances = {}
#     with open(filepath) as fh:
#         for line in fh:
#             line = line.strip()
#             if not line:
#                 continue
#             cols = line.split("\t")
#             if len(cols) < 3:
#                 continue
#             try:
#                 dist = float(cols[2])
#             except ValueError:
#                 continue
#             a = normalise_path(cols[0]).replace("plasmid_references", "binned_fragments")
#             b = normalise_path(cols[1]).replace("plasmid_references", "binned_fragments")
#             distances[(a, b)] = dist
#             distances[(b, a)] = dist  # store both orders
    
#     return distances
 
 
# def pair_label(paths):
#     ids = []
#     for p in paths:
#         m = re.search(r'(CAV\d+)', os.path.basename(p))
#         if m:
#             ids.append(m.group(1))
#     return " / ".join(ids) if ids else os.path.basename(paths[0])
 
 
# # ── plotting ──────────────────────────────────────────────────────────────────
 
# def plot_scatter(points):
#     x = np.array([p["mash_distance"] for p in points])
#     y = np.array([p["variance"] for p in points])
#     labels = [p["label"] for p in points]
 
#     fig, ax = plt.subplots(figsize=(9, 6))
 
#     scatter = ax.scatter(
#         x, y,
#         s=90, zorder=3,
#         c=y, cmap="viridis",
#         edgecolors="white", linewidths=0.6
#     )
 
#     # for label, xi, yi in zip(labels, x, y):
#     #     ax.annotate(
#     #         label, (xi, yi),
#     #         textcoords="offset points", xytext=(8, 4),
#     #         fontsize=8
#     #     )
 
#     cbar = fig.colorbar(scatter, ax=ax, pad=0.02)
#     cbar.set_label("Variance", fontsize=9)
 
#     # Spearman rank correlation + significance test
#     if len(x) >= 3:
#         r, p_value = stats.spearmanr(x, y)
#         r, p_value = float(r), float(p_value)
#         sig_label = "significant" if p_value < 0.05 else "not significant"
#         annotation = (
#             f"Spearman ρ = {r:.4f}\n"
#             f"p-value = {p_value:.4f}\n"
#             f"({sig_label} at α=0.05)\n"
#             f"n = {len(x)} pairs (subsample)"
#         )
#         # ax.text(
#         #     0.05, 0.95, annotation,
#         #     transform=ax.transAxes,
#         #     fontsize=9,
#         #     verticalalignment="top",
#         #     bbox=dict(boxstyle="round,pad=0.4", facecolor="lightyellow", edgecolor="gray", alpha=0.8)
#         # )
#         print(f"Spearman rank correlation (ρ): {r:.4f}")
#         print(f"p-value: {p_value:.4f} ({sig_label} at α=0.05)")
#         print(f"Note: based on subsample of {len(x)} pairs — interpret p-value with caution.")
#     elif len(x) == 2:
#         r, _ = stats.spearmanr(x, y)
#         # ax.text(
#         #     0.05, 0.95, f"Spearman ρ = {float(r):.4f}\n(need ≥3 pairs for p-value)",
#         #     transform=ax.transAxes,
#         #     fontsize=9, color="#ffa657",
#         #     verticalalignment="top",
#         #     bbox=dict(boxstyle="round,pad=0.4", facecolor="#1a1d27", edgecolor="#ffa657", alpha=0.8)
#         # )
#         print(f"Spearman ρ = {float(r):.4f} (need at least 3 pairs to compute p-value)")
#     else:
#         print("Not enough points to compute correlation (need at least 2).")
 
#     ax.set_xlabel("Mash distance between plasmid pair", fontsize=11)
#     ax.set_ylabel("Variance of distance distribution", fontsize=11)
#    # ax.set_title("Plasmid pair: mash distance vs bin variance", fontsize=13)
#     ax.grid(True, color="#30363d", linewidth=0.5, linestyle="--")
 
#     plt.tight_layout()
#     out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "binned_analysis.png")
#     plt.savefig(out_path, dpi=150, bbox_inches="tight")
#     print(f"Plot saved to: {out_path}")
#     plt.show()
 
 
# # ── main ──────────────────────────────────────────────────────────────────────
 
# def main():
#     if len(sys.argv) != 3:
#         print("Usage: python analyze_binned_fragments.py <binned_file> <mash_distance_file>")
#         sys.exit(1)
 
#     binned_file = sys.argv[1]
#     mash_file = sys.argv[2]
 
#     for path, label in [(binned_file, "binned file"), (mash_file, "mash distance file")]:
#         if not os.path.isfile(path):
#             print(f"Error: {label} not found: {path}")
#             sys.exit(1)
 
#     print(f"Parsing binned file: {binned_file}")
#     entries = parse_binned_file(binned_file)
#     print(f"Found {len(entries)} plasmid pair(s).\n")
 
#     print(f"Parsing mash distance file: {mash_file}")
#     mash_distances = parse_mash_file(mash_file)
#     print(f"Loaded {len(mash_distances) // 2} unique mash distance entries.\n")
 
#     points = []
#     for fasta_paths, numbers in entries:
#         #label = pair_label(fasta_paths)
#         print(label)
#         if len(fasta_paths) < 2:
#             print(f"  WARNING: pair [{label}] has fewer than 2 files. Skipping.")
#             continue
 
#         a = os.path.basename(fasta_paths[0])
#         b = os.path.basename(fasta_paths[1])

#         a = a.replace("_binned","")
#         b = b.replace("_binned","")
#         dist = mash_distances.get((a, b))
#         if dist is None:
#             print(f"  WARNING: no mash distance found for pair [{label}] ({a} vs {b}). Skipping.")
#             continue
 
#         variance = float(np.var(numbers))
#         print(f"  Pair [{label}]: mash distance = {dist:.6f}, variance = {variance:.4f}")
#         points.append({"label": label, "mash_distance": dist, "variance": variance})
 
#     if not points:
#         print("No valid plasmid pairs found. Cannot plot.")
#         sys.exit(1)
 
#     print(f"\nPlotting {len(points)} pairs...")
#     plot_scatter(points)
 
 
# if __name__ == "__main__":
#     main()

#### Here starts error vs plots

# """
# analyze_binned_fragments.py
 
# Takes two input files:
#   1. binned file: pairs of FASTA paths followed by a distribution of bin numbers
#   2. reference file: pairs of FASTA paths followed by a single integer
 
# For each plasmid pair matched by basename across the two files, computes:
#   error = |single_number - min(bin_distribution)|
 
# Then produces the same scatter plots as before but with error on the y-axis:
#   - error vs mash distance
#   - error vs total plasmid length (from a database folder)
#   - error vs length difference between two database folders
 
# Spearman rank correlation and p-value are shown on each plot.
 
# Usage:
#     python analyze_binned_fragments.py <binned_file> <reference_file> <mash_file> <db_folder_a> <db_folder_b>
 
# Example:
#     python analyze_binned_fragments.py binned.txt reference.txt mash.tsv db_a/ db_b/
# """
 
# import sys
# import os
# import re
# import numpy as np
# import matplotlib.pyplot as plt
# from scipy import stats
 
 
# # ── parsers ───────────────────────────────────────────────────────────────────
 
# def extract_basenames(paths):
#     """Return frozenset of basenames for matching across files."""
#     return frozenset(os.path.basename(p).replace("_binned","") for p in paths)
 
 
# def parse_fasta_paths_from_line(line):
#     """Split a line that may contain concatenated fasta paths."""
#     parts = re.split(r'(?=/)', line.strip())
#     return [p.strip() for p in parts if p.strip().endswith(".fasta")]
 
 
# def parse_binned_file(filepath):
#     """
#     Returns a dict: frozenset(basenames) -> list of bin numbers
#     """
#     entries = {}
#     current_paths = []
 
#     with open(filepath) as fh:
#         for line in fh:
#             line = line.strip()
#             if not line:
#                 continue
#             if ".fasta" in line:
#                 current_paths = parse_fasta_paths_from_line(line)
#             else:
#                 tokens = line.split()
#                 numbers = []
#                 for t in tokens:
#                     try:
#                         numbers.append(int(t))
#                     except ValueError:
#                         pass
#                 if numbers and current_paths:
#                     key = extract_basenames(current_paths)
#                     entries[key] = (current_paths, numbers)
#                     current_paths = []
#     return entries
 
 
# def parse_reference_file(filepath):
#     """
#     Returns a dict: frozenset(basenames) -> (paths, single_integer)
#     Basenames have 'plasmid_references' replaced with 'binned_fragments' for matching.
#     """
#     entries = {}
#     current_paths = []
 
#     with open(filepath) as fh:
#         for line in fh:
#             line = line.strip()
#             if not line:
#                 continue
#             if ".fasta" in line:
#                 current_paths = parse_fasta_paths_from_line(line)
#             else:
#                 tokens = line.split()
#                 if tokens and current_paths:
#                     try:
#                         value = int(tokens[0])
#                         # normalise basenames for matching
#                         normalised = [
#                             os.path.basename(p).replace("plasmid_references", "binned_fragments")
#                             for p in current_paths
#                         ]
#                         key = frozenset(normalised)
#                         entries[key] = (current_paths, value)
#                         current_paths = []
#                     except ValueError:
#                         pass
 
#     print(entries)
#     return entries
 
 
# def parse_mash_file(filepath):
#     """Returns dict: (basename_a, basename_b) -> mash_distance (both orders)."""
#     distances = {}
#     with open(filepath) as fh:
#         for line in fh:
#             cols = line.strip().split("\t")
#             if len(cols) < 3:
#                 continue
#             try:
#                 dist = float(cols[2])
#             except ValueError:
#                 continue
#             a = os.path.basename(cols[0]).replace("plasmid_references", "binned_fragments")
#             b = os.path.basename(cols[1]).replace("plasmid_references", "binned_fragments")
#             distances[(a, b)] = dist
#             distances[(b, a)] = dist
#     return distances
 
 
# def find_fasta_in_db(fasta_path, db_dir):
#     """
#     Walk db_dir to find a file whose name, after stripping _binned, matches
#     the basename of fasta_path (which already has _binned stripped).
#     Returns the full path of the actual file found.
#     """
#     target = os.path.basename(fasta_path).replace("_binned", "")
#     for root, dirs, files in os.walk(db_dir):
#         for fname in files:
#             if fname.replace("_binned", "") == target:
#                 return os.path.join(root, fname)
#     return None
 
 
# def total_sequence_length(filepath):
#     total = 0
#     with open(filepath) as fh:
#         for line in fh:
#             line = line.strip()
#             if line and not line.startswith(">"):
#                 total += len(line)
#     return total
 
 
# def avg_total_length_for_pair(fasta_paths, db_dir):
#     lengths = []
#     for fp in fasta_paths:
#         full_path = find_fasta_in_db(fp, db_dir)
#         if full_path is None:
#             return None
#         lengths.append(total_sequence_length(full_path))
#     return float(np.mean(lengths))
 
 
# def count_fasta_headers(filepath):
#     """Count lines starting with > in a FASTA file."""
#     count = 0
#     with open(filepath) as fh:
#         for line in fh:
#             if line.startswith(">"):
#                 count += 1
#     return count
 
 
# def avg_header_count_for_pair(fasta_paths, db_dir):
#     counts = []
#     for fp in fasta_paths:
#         full_path = find_fasta_in_db(fp, db_dir)
#         if full_path is None:
#             return None
#         counts.append(count_fasta_headers(full_path))
#     return float(np.mean(counts))
 
# def pair_label(paths):
#     ids = []
#     for p in paths:
#         m = re.search(r'(CAV\d+|[A-Z0-9]+(?:NIH|ECR)\w*)', os.path.basename(p))
#         if m:
#             ids.append(m.group(1))
#     if not ids:
#         ids = [os.path.basename(p)[:8] for p in paths[:2]]
#     return " / ".join(ids)
 
 
# # ── plotting ──────────────────────────────────────────────────────────────────
 
# def plot_scatter(x, y, labels, xlabel, title, outname):
#     fig, ax = plt.subplots(figsize=(9, 6))
 
#     scatter = ax.scatter(
#         x, y,
#         s=90, zorder=3,
#         c="black",
#         edgecolors="white", linewidths=0.6
#     )
 
#     # for label, xi, yi in zip(labels, x, y):
#     #     ax.annotate(label, (xi, yi),
#     #                 textcoords="offset points", xytext=(8, 4), fontsize=8)
 
#     # cbar = fig.colorbar(scatter, ax=ax, pad=0.02)
#     # cbar.set_label("Error", fontsize=9)
 
#     if len(x) >= 3:
#         r, p_value = stats.spearmanr(x, y)
#         r, p_value = float(r), float(p_value)
#         sig_label = "significant" if p_value < 0.05 else "not significant"
#         annotation = (
#             f"Spearman ρ = {r:.4f}\n"
#             f"p-value = {p_value:.4f}\n"
#             f"({sig_label} at α=0.05)\n"
#             f"n = {len(x)} pairs (subsample)"
#         )
#         # ax.text(0.05, 0.95, annotation, transform=ax.transAxes, fontsize=9,
#         #         verticalalignment="top",
#         #         bbox=dict(boxstyle="round,pad=0.4", facecolor="lightyellow", edgecolor="gray", alpha=0.8))
#         print(f"\n{title}")
#         print(f"  Spearman ρ: {r:.4f}, p-value: {p_value:.4f} ({sig_label} at α=0.05)")
#         print(f"  n = {len(x)} pairs")
 
#     ax.set_xlabel(xlabel, fontsize=11)
#     ax.set_ylabel("Error", fontsize=11)
#     ax.set_title(title, fontsize=13)
#     ax.grid(True, linestyle="--", alpha=0.5)
 
#     plt.tight_layout()
#     out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), outname)
#     plt.savefig(out_path, dpi=150, bbox_inches="tight")
#     print(f"  Plot saved to: {out_path}")
#     plt.show()
 
 
# # ── main ──────────────────────────────────────────────────────────────────────
 
# def main():
#     if len(sys.argv) != 6:
#         print("Usage: python analyze_binned_fragments.py <binned_file> <reference_file> <mash_file> <db_folder_a> <db_folder_b>")
#         sys.exit(1)
 
#     binned_file, ref_file, mash_file, folder_a, folder_b = sys.argv[1:]
 
#     for path, label in [(binned_file, "binned file"), (ref_file, "reference file"),
#                         (mash_file, "mash file"), (folder_a, "folder A"), (folder_b, "folder B")]:
#         if not os.path.exists(path):
#             print(f"Error: {label} not found: {path}")
#             sys.exit(1)
 
#     print("Parsing files...")
#     binned = parse_binned_file(binned_file)
#     reference = parse_reference_file(ref_file)
#     mash_distances = parse_mash_file(mash_file)
#     print(f"  Binned pairs: {len(binned)}")
#     print(f"  Reference pairs: {len(reference)}")
#     print(f"  Mash entries: {len(mash_distances) // 2}")
 
#     # Match pairs: reference keys are plain basenames, binned keys have _binned stripped
#     # so they should now match directly
#     matched = []
#     for ref_key, (ref_paths, single_number) in reference.items():
#         if ref_key not in binned:
#             print(f"  WARNING: no binned match for {pair_label(ref_paths)}. Skipping.")
#             continue
 
#         bin_paths, numbers = binned[ref_key]
#         error = abs(single_number - min(numbers))
#         label = pair_label(bin_paths)
 
#         # Mash lookup uses plain reference basenames
#         basenames = list(ref_key)
#         mash_dist = mash_distances.get((basenames[0], basenames[1]))
#         if mash_dist is None and len(basenames) >= 2:
#             mash_dist = mash_distances.get((basenames[1], basenames[0]))
 
#         avg_a = avg_total_length_for_pair(bin_paths, folder_a)
#         avg_b = avg_total_length_for_pair(bin_paths, folder_b)
#         length_diff = abs(avg_a - avg_b) if avg_a and avg_b else None
#         avg_total_len = avg_a if avg_a is not None else avg_b
#         avg_headers = avg_header_count_for_pair(bin_paths, folder_b)
 
#         matched.append({
#             "label": label,
#             "error": error,
#             "mash_dist": mash_dist,
#             "length_diff": length_diff,
#             "avg_total_length": avg_total_len,
#             "avg_header_count": avg_headers,
#         })
#         formatted = f"{mash_dist:.4f}" if mash_dist is not None else "N/A"
#         print(f"  [{label}]: ref={single_number}, min(bins)={min(numbers)}, error={error}, mash={formatted}")
 
#     print(f"\n{len(matched)} pairs matched successfully.")
 
#     # Plot 1: error vs mash distance
#     mash_points = [p for p in matched if p["mash_dist"] is not None]
#     if mash_points:
#         plot_scatter(
#             x=np.array([p["mash_dist"] for p in mash_points]),
#             y=np.array([p["error"] for p in mash_points]),
#             labels=[p["label"] for p in mash_points],
#             xlabel="Mash distance between plasmid pair",
#             title="",
#             outname="error_vs_mash.png"
#         )
 
#     # Plot 2: error vs length difference between folders
#     diff_points = [p for p in matched if p["length_diff"] is not None]
#     if diff_points:
#         plot_scatter(
#             x=np.array([p["length_diff"] for p in diff_points]),
#             y=np.array([p["error"] for p in diff_points]),
#             labels=[p["label"] for p in diff_points],
#             xlabel=f"Absolute difference in avg total length (bp)",
#             title="",
#             outname="error_vs_length_diff.png"
#         )
 
 
#     # Plot 3: error vs average total length
#     len_points = [p for p in matched if p["avg_total_length"] is not None]
#     if len_points:
#         plot_scatter(
#             x=np.array([p["avg_total_length"] for p in len_points]),
#             y=np.array([p["error"] for p in len_points]),
#             labels=[p["label"] for p in len_points],
#             xlabel="Average total plasmid length (bp)",
#             title="",
#             outname="error_vs_avg_length.png"
#         )
 
#     # Plot 4: error vs average number of fragments (>)
#     header_points = [p for p in matched if p["avg_header_count"] is not None]
#     if header_points:
#         plot_scatter(
#             x=np.array([p["avg_header_count"] for p in header_points]),
#             y=np.array([p["error"] for p in header_points]),
#             labels=[p["label"] for p in header_points],
#             xlabel="Average number of contigs",
#             title="",
#             outname="error_vs_avg_headers.png"
#         )
 
 
# if __name__ == "__main__":
#     main()
 

 ### pling vs mash of the pling dataset

# """
# Scatterplot: distance (file 1) vs. third-column ratio (file 2)
 
# File 1 format (TSV):  plasmid_1  plasmid_2  distance
# File 2 format (TSV):  plasmid_1  plasmid_2  col3  col4  col5
#   where col3 is a fraction like  1000/1000
# # """
 
# import argparse
# import sys
# import pandas as pd
# import matplotlib.pyplot as plt
# import os
# from scipy import stats
 
# # ── helpers ──────────────────────────────────────────────────────────────────
 
# def parse_fraction(value: str) -> float:
#     """Convert '950/1000' → 0.95.  Plain numbers pass through unchanged."""
#     value = str(value).strip()
#     if "/" in value:
#         num, den = value.split("/", 1)
#         return float(num) / float(den)
#     return float(value)
 
 
# def make_key(a: str, b: str) -> frozenset:
#     """Order-independent pair key so AP1–AP2 matches AP2–AP1."""
#     return frozenset([a.strip(), b.strip()])
 
 
# # ── loaders ──────────────────────────────────────────────────────────────────
 
# def load_distances(path: str) -> dict:
#     """Returns {frozenset(p1, p2): distance}."""
#     df = pd.read_csv(path, sep="\t", header=0,
#                      names=["plasmid_1", "plasmid_2", "distance"])
#     return {make_key(r.plasmid_1, r.plasmid_2): float(r.distance)
#             for _, r in df.iterrows()}
 
 
# def load_ratios(path: str) -> dict:
#     """Returns {frozenset(p1, p2): ratio} using the 3rd data column."""
#     df = pd.read_csv(path, sep="\t", header=None)
#     # Strip .fna suffixes if present so names match the distance file
#     df[0] = df[0].str.replace(r"\.fna$", "", regex=True)
#     df[1] = df[1].str.replace(r"\.fna$", "", regex=True)
#     print(df)
#     return {make_key(r[0], r[1]): float(r[2])
#             for _, r in df.iterrows() if r[2]<=0.15}
 
# def parse_mash_file(filepath):
#     """Returns dict: (basename_a, basename_b) -> mash_distance (both orders)."""
#     distances = []
#     with open(filepath) as fh:
#         for line in fh:
#             cols = line.strip().split("\t")
#             if len(cols) < 3:
#                 continue
#             try:
#                 dist = float(cols[2])
#             except ValueError:
#                 continue
#             a = os.path.basename(cols[0]).replace("plasmid_references", "binned_fragments")
#             b = os.path.basename(cols[1]).replace("plasmid_references", "binned_fragments")
#             distances.append = (frozenset(a,b),dist)
#             distances.append = (frozenset(a,b),dist)
#     return distances
# # ── main ─────────────────────────────────────────────────────────────────────
 
# def main():
#     parser = argparse.ArgumentParser(
#         description="Scatterplot: pairwise distance vs. alignment ratio"
#     )
#     parser.add_argument("distance_file",
#                         help="TSV with columns: plasmid_1, plasmid_2, distance")
#     parser.add_argument("ratio_file",
#                         help="TSV with columns: plasmid_1, plasmid_2, col3 (fraction), …")
#     parser.add_argument("-o", "--output", default="scatterplot.png",
#                         help="Output image path (default: scatterplot.png)")
#     parser.add_argument("--xlabel", default="Pairwise distance",
#                         help="X-axis label")
#     parser.add_argument("--ylabel", default="Alignment ratio (col 3)",
#                         help="Y-axis label")
#     parser.add_argument("--title", default="Distance vs. Alignment ratio",
#                         help="Plot title")
#     parser.add_argument("--alpha", type=float, default=0.7,
#                         help="Point transparency (0–1, default 0.7)")
#     parser.add_argument("--point-size", type=float, default=40,
#                         help="Marker size (default 40)")
#     args = parser.parse_args()
 
#     # Load both files
#     distances = load_distances(args.distance_file)
#     ratios    = load_ratios(args.ratio_file)
 
#     # Inner join on shared pairs
#     common_keys = set(distances) & set(ratios)
#     if not common_keys:
#         sys.exit(
#             "ERROR: No matching plasmid pairs found between the two files.\n"
#             "Check that names match (the script strips .fna suffixes automatically)."
#         )
 
#     x = [distances[k] for k in common_keys]
#     y = [ratios[k]    for k in common_keys]
 
#     print(f"Matched {len(common_keys)} pairs  "
#           f"({len(distances) - len(common_keys)} distance-only, "
#           f"{len(ratios) - len(common_keys)} ratio-only — skipped)")
 
#     # ── statistics ────────────────────────────────────────────────────────────
#     # Pearson (linear correlation)
#     pearson_r, pearson_p = stats.pearsonr(x, y)
 
#     # Spearman (rank-based, no normality assumption)
#     spearman_r, spearman_p = stats.spearmanr(x, y)
 
#     # Shapiro-Wilk normality test on both variables (only reliable for n < 5000)
#     sw_stat_x, sw_p_x = stats.shapiro(x) if len(x) <= 5000 else (float("nan"), float("nan"))
#     sw_stat_y, sw_p_y = stats.shapiro(y) if len(y) <= 5000 else (float("nan"), float("nan"))
 
#     def fmt_p(p):
#         if p < 0.001:
#             return "< 0.001"
#         return f"= {p:.3f}"
 
#     print("\n── Statistical results ──────────────────────────────")
#     print(f"  n pairs            : {len(x)}")
#     print(f"  Pearson  r         : {pearson_r:+.4f}   p {fmt_p(pearson_p)}")
#     print(f"  Spearman ρ         : {spearman_r:+.4f}   p {fmt_p(spearman_p)}")
#     print(f"  Shapiro-Wilk (x)   : W = {sw_stat_x:.4f},  p {fmt_p(sw_p_x)}")
#     print(f"  Shapiro-Wilk (y)   : W = {sw_stat_y:.4f},  p {fmt_p(sw_p_y)}")
#     print("─────────────────────────────────────────────────────\n")
 
#     # ── plot ──────────────────────────────────────────────────────────────────
#     #sns.set_theme(style="whitegrid", context="notebook")
 
#     fig, ax = plt.subplots(figsize=(8, 6))
 
#     ax.scatter(x, y,
#                s=args.point_size,
#                alpha=args.alpha,
#                edgecolors="white",
#                linewidths=0.4,
#                #color=sns.color_palette("mako", 1)[0]
#                )
 
#     # Regression line with 95 % CI
#     # sns.regplot(x=x, y=y, ax=ax, scatter=False,
#     #             line_kws={"color": "crimson", "linewidth": 1.5, "linestyle": "--"},
#     #             ci=95)
 
#     # Annotation box with both coefficients
#     annotation = (
#         f"Spearman ρ = {spearman_r:+.3f}  (p {fmt_p(spearman_p)})\n"
#         f"n = {len(x)}"
#     )
#     # ax.text(0.03, 0.97, annotation,
#     #         transform=ax.transAxes,
#     #         fontsize=9.5,
#     #         verticalalignment="top",
#     #         bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
#     #                   edgecolor="lightgray", alpha=0.9))
 
#     ax.set_xlabel("Mash Distance", fontsize=13)
#     ax.set_ylabel("Pling Distance", fontsize=13)
#     ax.set_title("", fontsize=15, fontweight="bold", pad=12)
 
#     plt.tight_layout()
#     plt.savefig(args.output, dpi=150)
#     print(f"Plot saved → {args.output}")
 
 
# if __name__ == "__main__":
#     main()