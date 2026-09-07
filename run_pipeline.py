"""
run_pipeline.py — entry point for starting the pipeline from a folder of FASTA files.

Takes a folder of plasmid FASTA files, builds pairs of them — either every
unique pair or a random subsample — and calls generateShortRead.py on each
pair in parallel batches, the same batching pattern groupCall.py's functions
use (a fixed-size chunk of subprocess.Popen calls, waited on before the next
chunk starts).

CLI:
    python run_pipeline.py --folder fastas --parallel 10
    python run_pipeline.py --folder fastas --parallel 10 --sample 50 --fragments 8 --implementation

Run from the Honour-Project directory (same assumption generateShortRead.py
and groupCall.py make about relative paths).
"""
import argparse
import itertools
import os
import random
import subprocess
import sys


def find_fasta_files(folder_path):
	exts = (".fasta", ".fna", ".fa")
	return sorted(f for f in os.listdir(folder_path) if f.lower().endswith(exts))


def build_pairs(folder, files, sample_size):
	"""Returns a list of (file1, file2) tuples, each formatted as
	'/<folder>/<name>' — the leading-slash form generateShortRead.py expects
	(it resolves them as '.' + file1). Uses every unique pair unless
	sample_size is given, in which case that many unique pairs are picked
	at random."""
	all_pairs = list(itertools.combinations(files, 2))
	if sample_size is None:
		chosen = all_pairs
	else:
		if sample_size > len(all_pairs):
			raise ValueError(
				f"--sample {sample_size} exceeds the number of unique pairs available ({len(all_pairs)})"
			)
		chosen = random.sample(all_pairs, sample_size)
	return [(f"/{folder}/{a}", f"/{folder}/{b}") for a, b in chosen]


def run_pairs(pairs, parallel, fragments, implementation):
	"""Runs generateShortRead.py on `pairs`, `parallel` at a time, waiting
	for each batch to finish before starting the next."""
	i = 0
	batch_id = 0
	while i < len(pairs):
		batch = pairs[i:i + parallel]
		commands = [
			f"python ./generateShortRead.py {file1} {file2} {batch_id} {index} {fragments} {implementation}"
			for index, (file1, file2) in enumerate(batch)
		]
		procs = [subprocess.Popen(cmd, shell=True) for cmd in commands]
		for p in procs:
			p.wait()
		print(f"The batch finished {i}")
		i += parallel
		batch_id += 1


def main():
	parser = argparse.ArgumentParser(
		description="Start the pipeline: pair up FASTA files from a folder and call generateShortRead.py on each pair."
	)
	parser.add_argument("--folder", required=True, help="Folder of FASTA files, relative to the current directory.")
	parser.add_argument("--parallel", type=int, default=10, help="Number of pairs to run concurrently per batch (default: 10).")
	parser.add_argument(
		"--sample", type=int, metavar="N", default=None,
		help="Randomly sample N unique pairs instead of using every pair in --folder."
	)
	parser.add_argument("--fragments", type=int, default=6, help="Used in the naming of files. Permits multiple instances of run_pipeline")
	parser.add_argument("--implementation", action="store_true", help="Use break_distance.my_implementation instead of pling.")
	args = parser.parse_args()

	folder = args.folder.strip("/")
	if not os.path.isdir(folder):
		sys.exit(f"Folder not found: {folder}")

	files = find_fasta_files(folder)
	if len(files) < 2:
		sys.exit(f"Need at least 2 FASTA files in {folder}, found {len(files)}")

	pairs = build_pairs(folder, files, args.sample)
	print(f"Running {len(pairs)} pair(s) from {folder}, {args.parallel} at a time")

	run_pairs(pairs, args.parallel, args.fragments, args.implementation)


if __name__ == "__main__":
	main()
