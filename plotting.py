"""
plotting.py — violin plots of per-pair distance distributions.

Reads an allDistances_v2_*.txt-style file (as written by
generateShortRead.writeTo()): pairs of lines where the first line names the
two plasmids compared and the second line is the space-separated distance
distribution sampled over fragment orderings. Renders one violin plot per
pair (or a red dot if the pair only has a single value) across paginated
PDF pages, annotated with each distribution's variance.

CLI: python plotting.py <input_file> [output_file] [pairs_per_page]
"""
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf_backend
import numpy as np
import sys
import os


def read(file):
    """Parses an allDistances_v2_*.txt-style file into
    {label_line: distance_array}, where label_line is the raw two-plasmid
    header line and distance_array the sampled distances beneath it."""
    data = {}
    distances = []
    comparisons = []
    groundTruth = []
    label = ""
    with open(file) as f:
        for line in f:
            if "/" in line:
                label = line
                #comparisons.append((line.split("/")[2].replace(".1.fna",""),line.split("/")[4].strip().replace(".1.fna","")))
                #groundTruth.append(float(subprocess.check_output(f'egrep "{line.split("/")[2].replace(".1.fna","")} {line.split("/")[4].strip().replace(".1.fna","")}" ./ground_Thruts.txt|cut -d " " -f3 ',shell=True,text=True).strip()))
            else:
                data[label] = np.fromstring(line.strip(),sep=" ")
                #distances.append(np.fromstring(line.strip(),sep=" "))
    return data

def load_distances(input_file):
    """Alternative parser to read(): same allDistances_v2_*.txt format, but
    keys the resulting dict with a "basename1\\nvs\\nbasename2" label
    (stripped of directory paths) instead of the raw header line, and
    splits/converts the distance line to floats up front. Not called by
    __main__ (which uses read()) — kept as a cleaner variant."""
    data = {}
    with open(input_file, "r") as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines) - 1:
        path_line = lines[i].strip()
        dist_line = lines[i+1].strip()
        i += 2

        files = path_line.split()
        if len(files) != 2:
            continue

        distances = [float(x) for x in dist_line.split() if x]
        label = f"{os.path.basename(files[0])}\nvs\n{os.path.basename(files[1])}"
        data[label] = distances

    return data

def plot_page(ax, labels, values):
    """Draws one page's worth of pairs onto `ax`: a violin plot for pairs
    with more than one sampled distance (annotated with variance), and a
    red scatter point for pairs with only a single value (variance=0)."""
    positions = range(len(labels))

    # Separate single and multi value entries
    single_pos  = [i for i, v in enumerate(values) if len(v) == 1]
    single_vals = [values[i][0] for i in single_pos]
    multi_pos   = [i for i, v in enumerate(values) if len(v) > 1]
    multi_vals  = [values[i] for i in multi_pos]

    # Violin for multi
    if multi_pos:
        parts = ax.violinplot(
            multi_vals,
            positions=multi_pos,
            showmedians=True,
            showextrema=True
        )
        for pc in parts["bodies"]:
            pc.set_alpha(0.7)

    # For each multi, annotate variance
    for i, pos in enumerate(multi_pos):
        var = np.var(values[pos])
        ax.annotate(
            f"var={var:.2f}",
            xy=(pos, max(values[pos])),
            fontsize=6,
            ha="center",
            va="bottom",
            color="blue"
        )

    # Scatter for single values
    if single_pos:
        ax.scatter(single_pos, single_vals, color="red", zorder=3, label="single value")
        for pos, val in zip(single_pos, single_vals):
            ax.annotate(
                f"var=0",
                xy=(pos, val),
                fontsize=6,
                ha="center",
                va="bottom",
                color="red"
            )

    ax.set_xticks(list(positions))
    ax.set_xticklabels(labels, fontsize=6, rotation=45, ha="right")
    ax.set_ylabel("Distance")
    ax.yaxis.grid(True, linestyle="--", alpha=0.5)
    ax.legend(fontsize=7)

def plot_distances(data, output_file, pairs_per_page=15):
    """Paginates `data` (as returned by read()/load_distances()) into
    groups of `pairs_per_page` and writes one page per group to a
    multi-page PDF via plot_page()."""
    labels = list(data.keys())
    values = list(data.values())

    # Split into pages
    pages = [
        (labels[i:i+pairs_per_page], values[i:i+pairs_per_page])
        for i in range(0, len(labels), pairs_per_page)
    ]

    with pdf_backend.PdfPages(output_file) as pdf:
        for page_num, (page_labels, page_values) in enumerate(pages):
            fig, ax = plt.subplots(figsize=(20, 8))
            ax.set_title(f"Plasmid pair distances — page {page_num + 1}/{len(pages)}")
            plot_page(ax, page_labels, page_values)
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
            print(f"[PAGE] {page_num + 1}/{len(pages)}")

    print(f"[SAVED] {output_file}")

if __name__ == "__main__":
    input_file      = sys.argv[1]
    output_file     = sys.argv[2] if len(sys.argv) > 2 else "distances_plot.pdf"
    pairs_per_page  = int(sys.argv[3]) if len(sys.argv) > 3 else 15

    data = read(input_file)
    plot_distances(data, output_file, pairs_per_page)