
"""
Takes two input files:
  1. binned file: pairs of FASTA paths followed by a distribution of bin numbers
  2. reference file: pairs of FASTA paths followed by a single integer
 
For each plasmid pair matched by basename across the two files, computes:
  error = |single_number - min(bin_distribution)|
 
Then produces the same scatter plots as before but with error on the y-axis:
  - error vs mash distance
  - error vs total plasmid length (from a database folder)
  - error vs length difference between two database folders
 
Spearman rank correlation and p-value are shown on each plot.
"""
 
import sys
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
 
 
# ── parsers ───────────────────────────────────────────────────────────────────
 
def extract_basenames(paths):
    """Return frozenset of basenames for matching across files."""
    return frozenset(os.path.basename(p).replace("_binned","") for p in paths)
 
 
def parse_fasta_paths_from_line(line):
    """Split a line that may contain concatenated fasta paths."""
    parts = re.split(r'(?=/)', line.strip())
    return [p.strip() for p in parts if p.strip().endswith(".fasta")]
 
 
def parse_binned_file(filepath):
    """
    Returns a dict: frozenset(basenames) -> list of bin numbers
    """
    entries = {}
    current_paths = []
 
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if ".fasta" in line:
                current_paths = parse_fasta_paths_from_line(line)
            else:
                tokens = line.split()
                numbers = []
                for t in tokens:
                    try:
                        numbers.append(int(t))
                    except ValueError:
                        pass
                if numbers and current_paths:
                    key = extract_basenames(current_paths)
                    entries[key] = (current_paths, numbers)
                    current_paths = []
    return entries
 
 
def parse_reference_file(filepath):
    """
    Returns a dict: frozenset(basenames) -> (paths, single_integer)
    Basenames have 'plasmid_references' replaced with 'binned_fragments' for matching.
    """
    entries = {}
    current_paths = []
 
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if ".fasta" in line:
                current_paths = parse_fasta_paths_from_line(line)
            else:
                tokens = line.split()
                if tokens and current_paths:
                    try:
                        value = int(tokens[0])
                        # normalise basenames for matching
                        normalised = [
                            os.path.basename(p).replace("plasmid_references", "binned_fragments")
                            for p in current_paths
                        ]
                        key = frozenset(normalised)
                        entries[key] = (current_paths, value)
                        current_paths = []
                    except ValueError:
                        pass
 
    print(entries)
    return entries
 
 
def parse_mash_file(filepath):
    """Returns dict: (basename_a, basename_b) -> mash_distance (both orders)."""
    distances = {}
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
            distances[(a, b)] = dist
            distances[(b, a)] = dist
    return distances
 
 
def find_fasta_in_db(fasta_path, db_dir):
    """
    Walk db_dir to find a file whose name, after stripping _binned, matches
    the basename of fasta_path (which already has _binned stripped).
    Returns the full path of the actual file found.
    """
    target = os.path.basename(fasta_path).replace("_binned", "")
    for root, dirs, files in os.walk(db_dir):
        for fname in files:
            if fname.replace("_binned", "") == target:
                return os.path.join(root, fname)
    return None
 
 
def total_sequence_length(filepath):
    total = 0
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith(">"):
                total += len(line)
    return total
 
 
def avg_total_length_for_pair(fasta_paths, db_dir):
    lengths = []
    for fp in fasta_paths:
        full_path = find_fasta_in_db(fp, db_dir)
        if full_path is None:
            return None
        lengths.append(total_sequence_length(full_path))
    return float(np.mean(lengths))
 
 
def count_fasta_headers(filepath):
    """Count lines starting with > in a FASTA file."""
    count = 0
    with open(filepath) as fh:
        for line in fh:
            if line.startswith(">"):
                count += 1
    return count
 
 
def avg_header_count_for_pair(fasta_paths, db_dir):
    counts = []
    for fp in fasta_paths:
        full_path = find_fasta_in_db(fp, db_dir)
        if full_path is None:
            return None
        counts.append(count_fasta_headers(full_path))
    return float(np.mean(counts))
 
def pair_label(paths):
    ids = []
    for p in paths:
        m = re.search(r'(CAV\d+|[A-Z0-9]+(?:NIH|ECR)\w*)', os.path.basename(p))
        if m:
            ids.append(m.group(1))
    if not ids:
        ids = [os.path.basename(p)[:8] for p in paths[:2]]
    return " / ".join(ids)
 
 
# ── plotting ──────────────────────────────────────────────────────────────────
 
def plot_scatter(x, y, labels, xlabel, title, outname):
    fig, ax = plt.subplots(figsize=(9, 6))
 
    scatter = ax.scatter(
        x, y,
        s=90, zorder=3,
        c="black",
        edgecolors="white", linewidths=0.6
    )
 
    # for label, xi, yi in zip(labels, x, y):
    #     ax.annotate(label, (xi, yi),
    #                 textcoords="offset points", xytext=(8, 4), fontsize=8)
 
    # cbar = fig.colorbar(scatter, ax=ax, pad=0.02)
    # cbar.set_label("Error", fontsize=9)
 
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
        # ax.text(0.05, 0.95, annotation, transform=ax.transAxes, fontsize=9,
        #         verticalalignment="top",
        #         bbox=dict(boxstyle="round,pad=0.4", facecolor="lightyellow", edgecolor="gray", alpha=0.8))
        print(f"\n{title}")
        print(f"  Spearman ρ: {r:.4f}, p-value: {p_value:.4f} ({sig_label} at α=0.05)")
        print(f"  n = {len(x)} pairs")
 
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel("Error", fontsize=11)
    ax.set_title(title, fontsize=13)
    ax.grid(True, linestyle="--", alpha=0.5)
 
    plt.tight_layout()
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), outname)
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"  Plot saved to: {out_path}")
    plt.show()
 
 
# ── main ──────────────────────────────────────────────────────────────────────
 
def main():
    if len(sys.argv) != 6:
        print("Usage: python analyze_binned_fragments.py <binned_file> <reference_file> <mash_file> <db_folder_a> <db_folder_b>")
        sys.exit(1)
 
    binned_file, ref_file, mash_file, folder_a, folder_b = sys.argv[1:]
 
    for path, label in [(binned_file, "binned file"), (ref_file, "reference file"),
                        (mash_file, "mash file"), (folder_a, "folder A"), (folder_b, "folder B")]:
        if not os.path.exists(path):
            print(f"Error: {label} not found: {path}")
            sys.exit(1)
 
    print("Parsing files...")
    binned = parse_binned_file(binned_file)
    reference = parse_reference_file(ref_file)
    mash_distances = parse_mash_file(mash_file)
    print(f"  Binned pairs: {len(binned)}")
    print(f"  Reference pairs: {len(reference)}")
    print(f"  Mash entries: {len(mash_distances) // 2}")
 
    # Match pairs: reference keys are plain basenames, binned keys have _binned stripped
    # so they should now match directly
    matched = []
    for ref_key, (ref_paths, single_number) in reference.items():
        if ref_key not in binned:
            print(f"  WARNING: no binned match for {pair_label(ref_paths)}. Skipping.")
            continue
 
        bin_paths, numbers = binned[ref_key]
        error = abs(single_number - min(numbers))
        label = pair_label(bin_paths)
 
        # Mash lookup uses plain reference basenames
        basenames = list(ref_key)
        mash_dist = mash_distances.get((basenames[0], basenames[1]))
        if mash_dist is None and len(basenames) >= 2:
            mash_dist = mash_distances.get((basenames[1], basenames[0]))
 
        avg_a = avg_total_length_for_pair(bin_paths, folder_a)
        avg_b = avg_total_length_for_pair(bin_paths, folder_b)
        length_diff = abs(avg_a - avg_b) if avg_a and avg_b else None
        avg_total_len = avg_a if avg_a is not None else avg_b
        avg_headers = avg_header_count_for_pair(bin_paths, folder_b)
 
        matched.append({
            "label": label,
            "error": error,
            "mash_dist": mash_dist,
            "length_diff": length_diff,
            "avg_total_length": avg_total_len,
            "avg_header_count": avg_headers,
        })
        formatted = f"{mash_dist:.4f}" if mash_dist is not None else "N/A"
        print(f"  [{label}]: ref={single_number}, min(bins)={min(numbers)}, error={error}, mash={formatted}")
 
    print(f"\n{len(matched)} pairs matched successfully.")
 
    # Plot 1: error vs mash distance
    mash_points = [p for p in matched if p["mash_dist"] is not None]
    if mash_points:
        plot_scatter(
            x=np.array([p["mash_dist"] for p in mash_points]),
            y=np.array([p["error"] for p in mash_points]),
            labels=[p["label"] for p in mash_points],
            xlabel="Mash distance between plasmid pair",
            title="",
            outname="error_vs_mash.png"
        )
 
    # Plot 2: error vs length difference between folders
    diff_points = [p for p in matched if p["length_diff"] is not None]
    if diff_points:
        plot_scatter(
            x=np.array([p["length_diff"] for p in diff_points]),
            y=np.array([p["error"] for p in diff_points]),
            labels=[p["label"] for p in diff_points],
            xlabel=f"Absolute difference in avg total length (bp)",
            title="",
            outname="error_vs_length_diff.png"
        )
 
 
    # Plot 3: error vs average total length
    len_points = [p for p in matched if p["avg_total_length"] is not None]
    if len_points:
        plot_scatter(
            x=np.array([p["avg_total_length"] for p in len_points]),
            y=np.array([p["error"] for p in len_points]),
            labels=[p["label"] for p in len_points],
            xlabel="Average total plasmid length (bp)",
            title="",
            outname="error_vs_avg_length.png"
        )
 
    # Plot 4: error vs average number of fragments (>)
    header_points = [p for p in matched if p["avg_header_count"] is not None]
    if header_points:
        plot_scatter(
            x=np.array([p["avg_header_count"] for p in header_points]),
            y=np.array([p["error"] for p in header_points]),
            labels=[p["label"] for p in header_points],
            xlabel="Average number of contigs",
            title="",
            outname="error_vs_avg_headers.png"
        )
 
 
if __name__ == "__main__":
    main()