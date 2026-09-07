"""
generateShortRead.py — core distance-sampling driver.

Given two already-fragmented plasmid assemblies (each fragment is a separate
FASTA record, standing in for contigs assembled from short reads), this
script reconstructs both plasmids under different orderings of their
fragments and measures the rearrangement distance between them for each
ordering. Sampling over many fragment orderings gives a *distribution* of
distances for one plasmid pair rather than a single number, because the
true contig order is unknown from short-read assembly alone.

Two distance backends are supported (selected by the `implementation` CLI
flag):
  - pling (external tool, called via subprocess) — containment + DCJ
    distance, invoked through the `pling` CLI.
  - `break_distance.my_implementation` — this project's own nucmer-based
    synteny-block + 2-break distance implementation.

Entry point: run as `python generateShortRead.py <file1> <file2> <iteration>
<parallel_id> <nrOfFragments> <implementation>`. Intended to be invoked in
batches by groupCall.py rather than run standalone.
"""
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import sys
import random as rdm
from itertools import permutations
from break_distance import my_implementation

def writeTo(resultFile,distances,file1,file2):
	"""Append one result block to resultFile: a header line with the two
	input filenames concatenated, followed by the sampled distances
	space-separated on the next line."""
	resultFile.write(file1+file2+"\n")
	for i in distances:
		resultFile.write(str(i)+" ")
	resultFile.write("\n")

def generateShortOld(nr,file,nrOfFragments) :
	"""Legacy fragmenter: splits a single-record FASTA into nrOfFragments
	roughly-equal chunks by raw line count, rather than reading pre-existing
	fragment records. Superseded by generateShort(); kept for reference."""
	i = 0
	with open(file) as f:
		val = []
		val.append(f.readline().strip())
		s = ""
		for x in f:
			s+=x.strip()
			i+=1
			if i == int(nr/nrOfFragments)+1:
				val.append(s)
				s=""
				i=0
		val.append(s)
		return val

def generateShort(nr, file, nrOfFragments):
	"""
	Reads an already-fragmented assembly FASTA file.
	Each fragment is a separate sequence entry (with its own header).
	Returns a list where index 0 is unused (kept for structural compatibility),
	and indices 1..N are the concatenated sequence lines of each fragment.
	"""
	fragments = []
	current_seq = []

	with open(file) as f:
		for line in f:
			line = line.strip()
			if line.startswith(">"):
				if current_seq:
					fragments.append("".join(current_seq))
					current_seq = []
				# Skip storing headers — construct() uses val[0] for the single
				# output header, so individual fragment headers are discarded
			else:
				current_seq.append(line)
		if current_seq:
			fragments.append("".join(current_seq))

	# val[0] is the FASTA header written by construct() into the output file.
	# We derive it from the filename to keep a meaningful label.
	header = ">" + file.split("/")[-1].replace(".fasta", "")
	return [header] + fragments

def try_something(val):
	"""Debug/scratch helper: writes each fragment of val on its own as a
	single-record FASTA and runs pling against it individually. Not used by
	the main pipeline (initial/initial_my_implementation)."""
	distances = []
	for i in range(1,len(val)):
		newS = ""
		newS+=val[i]
		g = open(f"./py_out/try_something.fasta","w")
		g.write(val[0]+'\n'+newS+'\n')
		g.close()

		subprocess.call(f"pling $PWD/py_out/input_test_something.txt output_dir_something_{i} align --containment_distance 1 --dcj 10",shell=True)
		distances.append(int(distance(f"output_dir_something_{i}")))
		subprocess.call(f"rm -r ./output_dir_something_{i}", shell=True)
	return distances

def distance(s):
	"""Reads pling's output TSV from directory s (all_plasmids_distances.tsv)
	and returns the distance value from the second row, third column.
	Returns -1 if pling produced no result row (e.g. it failed to align)."""
	#good reminder I quess, not exactly necessary but in case pling does not calculate the distance we add -1
	f = open(s+"/all_plasmids_distances.tsv")
	f.readline()
	values = f.readline()
	if values != "":
		return values.split('\t')[2].strip()
	else:
		return -1

def construct(indexes, iteration, parralel, fragments, suffix):
	"""
	Reconstructs a FASTA file from the given fragment list using the provided index permutation.
	suffix distinguishes file1 vs file2 outputs (e.g. 'a' and 'b').
	"""
	s = ""
	for i in indexes:
		s += fragments[i + 1]
	i = 0
	newS = ""
	while i < len(s):
		newS += s[i:i+70] + '\n'
		i += 70
	g = open(f"./py_out/try{iteration}_{parralel}_{suffix}.fasta", "w")
	g.write(fragments[0] + '\n' + newS + '\n')
	g.close()
	
def generate(k,size,iteration,used,elements,parralel):
	"""Exhaustively generates all permutations and calculates pling distances."""
	if k == size:
		global distances
		global var
		var += 1
		construct(elements, iteration, parralel, val1, "a")
		construct(elements, iteration, parralel, val2, "b")
		subprocess.call(
			f"pling $PWD/py_out/input{iteration}_{parralel}.txt output_dir_{var}_{parralel}_{size} align --containment_distance 1 --dcj 10",
			shell=True
		)
		distances.append(int(distance(f"output_dir_{var}_{parralel}_{size}")))
		subprocess.call(f"rm -r ./output_dir_{var}_{parralel}_{size}", shell=True)
	else:
		for i in range(size):
			if used[i] == False:
				elements[k] = i
				used[i] = True
				generate(k + 1, size, iteration, parralel, used, elements)
				used[i] = False

def factorial(size1,size2,target):
	"""Returns True if size1! * size2! <= target, False otherwise — i.e.
	whether the full space of (permutation-of-fragments-1, permutation-of-
	fragments-2) pairs is small enough to enumerate exhaustively rather than
	random-sample. Bails out early via the running-product check to avoid
	overflow/slowness on large inputs. Despite the name, this does not
	return a factorial value."""
	prod = 1
	for i in range(1,size1+1):
		if prod>target:
			return False
		prod*=i
	if prod>target:
			return False
	for i in range(1,size2+1):
		if prod>target:
			return False
		prod*=i
	if prod>target:
			return False
	return True

def reverse(n,size):
	"""Decodes integer n (0-indexed) into the n-th permutation of range(size)
	in Lehmer-code order — the standard "factorial number system" trick for
	turning a random integer directly into a permutation without generating
	the whole permutation list. Used only by the legacy subsample_Old();
	subsample() uses random.sample/itertools.permutations instead."""
	#this needs comments so here I go
	elems = [0 for i in range(size)]
	used = [False for i in range(size+1)]
	k = 0
	nextSmallest = 1
	#the main idea is to best fit n
	while (n != 0):
		#the if finds out if we put the next smallest or we have another number
		if factorial(size - k - 1) <= n:
			#we get the next nr
			elems[k] = get(int(n / factorial(size - k - 1)),used)
			
			#if there is no next nr
			if elems[k] == -1:
				break
			
			#tricky part here is that we want the left unused parts of n
			n = n % factorial(size - k - 1)
			
			used[elems[k]] = True
			k+=1
		else :
			#put the smallest and then find the next smallest
			elems[k] = nextSmallest
			k+=1
			used[nextSmallest] = True
			nextSmallest = nextSmall(used)
	if n > 0:
		return None
	#there is the chance that while we have no n left, we still have spaces in the array, so we just fill it
	elems = fill(elems, k,used)
	empty(used)
	return elems

# --- helpers for reverse() (Lehmer-code decoding) ---
def empty(used):
	"""Resets the `used` marker array back to all-False, in place."""
	for i in range(1,len(used)):
		used[i] = False
	return used
def fill(elems, k,used):
	"""Appends the still-unused values, in ascending order, into elems
	starting at index k — used to pad out a partially-decoded permutation."""
	for i in range(1,len(used)):
		if used[i]==False:
				elems[k] = i
				k+=1
	return elems
def get(n,used):
	"""Returns the (n+1)-th value not yet marked used, or -1 if there is
	no such value."""
	k = 0;
	for i in range(1,len(used)):
		if used[i]==False:
			k+=1
		if k == n + 1:
			return i
	return -1
def nextSmall(used):
	"""Returns the smallest value not yet marked used."""
	for i in range(1,len(used)):
		if used[i]==False:
				return i;
	return len(used) - 1;

def subsample_Old(size,file1,file2, iteration,parralel,implementation=False):
	"""Legacy sampler: draws 120 random permutation-pair indices via
	reverse() and runs pling/my_implementation on each reconstruction.
	Superseded by subsample(), which switches between exhaustive
	enumeration and random sampling depending on the size of the
	permutation space."""
	distances_sample = []
	fact = factorial(size)
	for i in range(120):
		index = rdm.randint(0,fact-1)
		elements = reverse(index,size)
		elements = [x-1 for x in elements]
		construct(elements,iteration,parralel)

		if implementation:
			distances_sample.append(my_implementation(file2,f"./py_out/try{iteration}_{parralel}.fasta"))
		else:
		# I put parralel so different threads do not fight between themselves
		# I put sizes so that same parralel counts for different fragment sized do not fight amongst themselves
		# I put i so that same parralel calls do not fight against themselves
			subprocess.run(f"pling $PWD/py_out/input{iteration}_{parralel}.txt output_dir_{i}_{parralel}_{size} align --containment_distance 1 --dcj 10",shell=True,check=True)
			distances_sample.append(int(distance(f"output_dir_{i}_{parralel}_{size}")))
			subprocess.call(f"rm -r ./output_dir_{i}_{parralel}_{size}", shell=True)

	print(len(distances_sample))
	w = open(f"allDistances_v2_{size}_{implementation}_{parralel}.txt","a")
	writeTo(w,distances_sample,file1,file2)
	w.close()

def subsample(size, file1, file2, iteration, parralel, implementation=False):
	"""Samples up to `target` (120) (permutation-of-val1, permutation-of-val2)
	pairs, reconstructs both plasmids for each pair, and measures the
	distance between them (pling or my_implementation depending on the
	`implementation` flag). If the full permutation space (size1! * size2!)
	is small enough, all pairs are enumerated exhaustively; otherwise pairs
	are drawn uniformly at random without replacement. Appends the
	resulting distance distribution to allDistances_v2_<size>_<implementation>_<parralel>.txt.
	Relies on globals val1/val2 set by initial()."""

	distances_sample = []
	size1 = len(val1) - 1
	size2 = len(val2) - 1
	indices1 = list(range(size1))
	indices2 = list(range(size2))
	target = 120

	if factorial(size1,size2,target):
		# Exhaustively enumerate all pairs
		all_pairs = [(p1, p2) for p1 in permutations(indices1) for p2 in permutations(indices2)]
		sampled_pairs = all_pairs
	else:
		# Random sampling for large spaces
		seen = set()
		sampled_pairs = []
		while len(sampled_pairs) < target:
			perm1 = tuple(rdm.sample(indices1, size1))
			perm2 = tuple(rdm.sample(indices2, size2))
			pair = (perm1, perm2)
			if pair not in seen:
				seen.add(pair)
				sampled_pairs.append(pair)

	for perm1, perm2 in sampled_pairs:
		construct(list(perm1), iteration, parralel, val1, "a")
		construct(list(perm2), iteration, parralel, val2, "b")

		if implementation:
			distances_sample.append(my_implementation(
				f"./py_out/try{iteration}_{parralel}_b.fasta",
				f"./py_out/try{iteration}_{parralel}_a.fasta"
			))
		else:
			i = len(distances_sample)
			subprocess.run(
				f"pling $PWD/py_out/input{iteration}_{parralel}.txt output_dir_{i}_{parralel}_{size} align --containment_distance 1 --dcj 10",
				shell=True, check=True
			)
			distances_sample.append(int(distance(f"output_dir_{i}_{parralel}_{size}")))
			subprocess.call(f"rm -r ./output_dir_{i}_{parralel}_{size}", shell=True)

	print(len(distances_sample))
	w = open(f"allDistances_v2_{size}_{implementation}_{parralel}.txt", "a")
	writeTo(w, distances_sample, file1, file2)
	w.close()


def initial(file1,file2,iteration,parralel,nrOfFragments,implementation):
	"""Pipeline entry point for one plasmid pair. Loads both fragmented
	FASTA files into globals val1/val2, writes the pling input-pair file
	for this iteration/parallel slot, then calls subsample() to sample
	fragment orderings and record the resulting distance distribution."""

	nr1 = subprocess.check_output(f"wc -l .{file1}", shell=True, text=True).strip().split(" ")[0]
	nr2 = subprocess.check_output(f"wc -l .{file2}", shell=True, text=True).strip().split(" ")[0]

	path = subprocess.check_output("echo $PWD", shell=True, text=True).strip()

	global val1, val2
	val1 = generateShort(int(nr1), "." + file1, nrOfFragments)
	val2 = generateShort(int(nr2), "." + file2, nrOfFragments)

	global distances
	distances = []

	g = open(f"./py_out/input{iteration}_{parralel}.txt", "w")
	g.write(path + f"/py_out/try{iteration}_{parralel}_a.fasta\n")
	g.write(path + f"/py_out/try{iteration}_{parralel}_b.fasta")
	g.close()

	subsample(nrOfFragments, file1, file2, iteration=iteration, parralel=parralel, implementation=implementation)
	
def initial_my_implementation(file1,file2,iteration,parralel,nrOfFragments):
	"""Variant of initial() for the single-plasmid case used with
	break_distance.my_implementation: fragments one plasmid (val) rather
	than a pair, then samples orderings via subsample(..., implementation=True)."""
	#The number of lines of the fasta
	nr = subprocess.check_output(f"wc -l .{file1}",shell=True,text=True)
	nr = nr.strip().split(" ")[0]

	path = subprocess.check_output("echo $PWD",shell=True,text=True).strip()

	#this generates a list of plasmid fragments
	global val
	val = generateShort(int(nr),"."+file1,nrOfFragments)

	global distances
	distances = []

	subsample(nrOfFragments,file1,file2,iteration,parralel,implementation=True)
	

if __name__=="__main__":
	# CLI: python generateShortRead.py <file1> <file2> <iteration> <parralel> <nrOfFragments> <implementation:true|false>
	file1 = sys.argv[1]
	file2 = sys.argv[2]
	iteration = int(sys.argv[3])
	parralel = int(sys.argv[4])
	nrOfFragments = int(sys.argv[5])
	implementation = sys.argv[6].lower() == "true"
	if implementation:
		initial_my_implementation(file1,file2,iteration,parralel,nrOfFragments)
	else:
		initial(file1,file2,iteration,parralel,nrOfFragments,implementation)


	