"""
bin_fragments.py — assigns real short-read contigs to their reference plasmid.

Takes a fragmentedDatabase/ assembly (short-read contigs for one sample,
identified by an SRR/DRR accession) and, using an SRR-id,
aligns those contigs against every plasmid_references/ sequence
belonging to that organism with nucmer. Each contig is then "binned" into
one or more output FASTA files under binned_fragments/ — one per reference
plasmid it aligns to with >85% coverage — via bin_fragments(). This is what
builds the binned_fragments/ dataset used elsewhere in the pipeline as the
real-short-read test case (see groupCall.py's data-source notes).

CLI: python bin_fragments.py <mapping_file> <plasmid_dir> <output_dir>
  Iterates every file in ./fragmentedDatabase and bins it.
"""
import subprocess
import os
import sys
import re
import collections

def load_mapping(mapping_file):
    """Reads a two-column whitespace-separated file into a dict:
    {SRR/DRR accession: organism name}."""
    mapping = {}
    with open(mapping_file, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                mapping[parts[0]] = parts[1]
    return mapping

def run_nucmer(reference, query, output_prefix):
    """Runs the standard nucmer -> delta-filter -> show-coords chain,
    aligning `query` against `reference` and writing
    {output_prefix}.{delta,filtered.delta,coords}."""
    subprocess.run(
        f"nucmer --maxmatch --prefix={output_prefix} {reference} {query}",
        shell=True, check=True
    )
    subprocess.run(
        f"delta-filter -q {output_prefix}.delta > {output_prefix}.filtered.delta",
        shell=True, check=True
    )
    subprocess.run(
        f"show-coords -rcl {output_prefix}.filtered.delta > {output_prefix}.coords",
        shell=True, check=True
    )

def merge_intervals(intervals):
    """
    Given a list of (start, end) tuples, merge overlapping/adjacent intervals
    and return the total covered bases.
    """
    if not intervals:
        return 0
    sorted_ivs = sorted(intervals, key=lambda x: x[0])
    merged = [sorted_ivs[0]]
    for start, end in sorted_ivs[1:]:
        if start <= merged[-1][1]:          # overlapping or adjacent
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return sum(end - start+1 for start, end in merged)

def parse_coords(coords_file,org_name):
    """
    Parses a show-coords file (reference=combined plasmid FASTA, query=the
    fragmented-contig FASTA) into per-fragment, per-plasmid coverage
    fractions. Returns scores[plasmid][fragment] = (merged aligned bases of
    that fragment against that plasmid) / (plasmid length) — i.e. how much
    of the *reference plasmid* this one fragment covers. bin_fragments()
    then bins a fragment into every plasmid where this exceeds 0.85.
    Overlapping/adjacent alignment intervals for the same (plasmid,
    fragment) pair are merged via merge_intervals() before computing
    coverage, so fragments aligning in multiple pieces aren't over-counted.
    """
    scores = collections.defaultdict(lambda: collections.defaultdict(float))  # {query_fragment: {plasmid_file: total_aligned_bases}}
    intervals = collections.defaultdict(lambda: collections.defaultdict(list))
    lengths = collections.defaultdict(lambda: collections.defaultdict(float))
    with open(coords_file, "r") as f:
        # Skip the header lines
        for _ in range(5):
            f.readline()
        for line in f:
            parts = line.strip().split("|")
            
            if len(parts) < 5:
                continue

            start_end = parts[1].strip()
            s = float(re.split(" +",start_end)[0])
            e = float(re.split(" +",start_end)[1])
            len2 = parts[2].strip()
            len2 = float(re.split(" +",len2)[1])
            qup = parts[6].strip()
            query_fragment = qup.split("\t")[0]
            plasmid        = qup.split("\t")[1]
            og_length = parts[4].strip()
            lenq = float(re.split(" +",og_length)[1])

            s,e = min(s, e), max(s, e)

            intervals[plasmid][query_fragment].append((s, e))
            # be aware that query fragment here is actually the reference and plasmid is the query
            if plasmid not in lengths:
                lengths[plasmid] = {}
            lengths[plasmid] = lenq

    for plasmid,fragments in intervals.items():
        for fragment,inter in fragments.items():
            scores[plasmid][fragment] = merge_intervals(inter)/lengths[plasmid]


    #For each fragment pick the plasmid with the most aligned bases
    # best = {}
    # for fragment, plasmid_counts in matches.items():
    #     best[fragment] = max(plasmid_counts, key=plasmid_counts.get)
    return scores

def read_fasta(fasta_file):
    """Returns a dict of {header: sequence}"""
    sequences = {}
    current_header = None
    current_seq = []
    with open(fasta_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header is not None:
                    sequences[current_header] = "".join(current_seq)
                current_header = line[1:]  # strip >
                current_seq = []
            else:
                current_seq.append(line)
        if current_header is not None:
            sequences[current_header] = "".join(current_seq)
    return sequences

def write_fasta(sequences, output_file):
    with open(output_file, "a") as f:
        for header, seq in sequences.items():
            f.write(f">{header}\n")
            i = 0
            while i < len(seq):
                f.write(seq[i:i+70] + "\n")
                i += 70

def bin_fragments(fasta_file, mapping_file, plasmid_dir, output_dir):
    """Bins every contig in `fasta_file` (short-read assembly for one
    sample) into whichever plasmid_dir reference sequence(s) it aligns to
    with >85% reference coverage (see parse_coords()). A fragment can be
    written into more than one plasmid's binned FASTA if it aligns well to
    several references (e.g. shared/repeated regions). Writes
    {plasmid_name}_binned.fasta files under output_dir. The organism name
    is derived from the SRR/DRR id embedded in fasta_file's filename via
    mapping_file (load_mapping()), and only that organism's reference
    plasmids in plasmid_dir are considered."""
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs("./nucmer_tmp", exist_ok=True)

    # Extract ID from filename e.g. pspades_DRR016448.fasta -> DRR016448
    basename = os.path.splitext(os.path.basename(fasta_file))[0]
    srr_id   = basename.split("_", 1)[1]

    # Map SRR/DRR to organism name e.g. AVNIH1
    mapping  = load_mapping(mapping_file)
    org_name = mapping.get(srr_id)
    if not org_name:
        print(f"[ERROR] No mapping found for {srr_id}")
        return

    # Find all plasmid files for this organism e.g. AVNIH1_1.fasta, AVNIH1_2.fasta ...
    plasmid_files = [
        os.path.join(plasmid_dir, f)
        for f in os.listdir(plasmid_dir)
        if f.startswith(org_name) and f.endswith(".fasta")
    ]

    if not plasmid_files:
        print(f"[ERROR] No plasmid files found for {org_name} in {plasmid_dir}")
        return

    print(f"[INFO] Found {len(plasmid_files)} plasmids for {org_name}")

    # Read all fragments
    fragments = read_fasta(fasta_file)

    assigned = {}  # {fragment_header: plasmid_file}

    combined_plasmids = f"./nucmer_tmp/{srr_id}_all_plasmids.fasta"
    with open(combined_plasmids, "w") as out:
        for pf in plasmid_files:
            with open(pf) as f:
                out.write(f.read())

    output_prefix = f"./nucmer_tmp/{srr_id}_all_vs_fragments"

    run_nucmer(combined_plasmids, fasta_file, output_prefix)

    coords_file = f"{output_prefix}.coords"
    scores = parse_coords(coords_file,org_name)
    assigned = {}

    for fragment_header, fragment_seq in fragments.items():
        if fragment_header not in scores:
            print(f"[UNASSIGNED] {fragment_header}")
            continue

        ref_scores = scores[fragment_header]
        #best_ref = max(ref_scores, key=ref_scores.get)

        for pl_header, score in ref_scores.items():
            if fragment_header not in assigned:
                assigned[fragment_header] = []
            if score>0.85:
                assigned[fragment_header].append(pl_header)

        for best_ref in assigned[fragment_header]:
            plasmid_name = os.path.splitext(os.path.basename(best_ref))[0]
            out_file = os.path.join(output_dir, f"{plasmid_name}_binned.fasta")

            write_fasta({fragment_header: fragment_seq}, out_file)

            print(f"[BINNED] {fragment_header} -> {plasmid_name}")
    
    subprocess.call("rm -r ./nucmer_tmp", shell=True)


def get_fragment_coverage(coords_file):
    """
    Sums up the query aligned bases (LEN 2) across all alignments,
    giving the total coverage of the fragment.
    """
    total = 0
    with open(coords_file, "r") as f:
        for _ in range(5):
            f.readline()
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("|")
            # if len(parts) < 2:
            #     continue
            # numeric_cols = parts[0].strip().split("|")
            # len_cols     = numeric_cols[2].strip().split()
            total       += int(re.split(" +",parts[2])[2])  # LEN 2 = query alignment length
    return total

if __name__ == "__main__":
    # Bins every assembly in ./fragmentedDatabase against plasmid_dir,
    # writing binned_fragments/-style output under output_dir.
    #fasta_file   = sys.argv[1]   # e.g. pspades_DRR016448.fasta
    mapping_file = sys.argv[1]   # all_genomes_files.txt
    plasmid_dir  = sys.argv[2]   # directory with CAV1311_1.fasta etc.
    output_dir   = sys.argv[3]   # where to write the binned fragments

    files = subprocess.check_output("ls -1 ./fragmentedDatabase",shell=True,text=True)
    files = files.split("\n")
    for file in files:
        bin_fragments("./fragmentedDatabase/"+file, mapping_file, plasmid_dir, output_dir)