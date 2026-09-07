"""
run_pling.py — baseline (no fragment-order sampling) pling distance.

Where generateShortRead.py samples many fragment orderings of a plasmid
pair to get a distance *distribution*, this module just runs pling once per
call and repeats that fixed call 5 times to check pling's own run-to-run
variability. Used by groupCall.run_base_pling_on_syntethic() as a baseline
to compare against the fragment-order-sampling results.
"""
import subprocess
def writeTo(resultFile,distances,file1,file2):
	"""Append one result block to resultFile: a header line with the two
	input filenames concatenated, followed by the sampled distances
	space-separated on the next line."""
	resultFile.write(file1+file2+"\n")
	for i in distances:
		resultFile.write(str(i)+" ")
	resultFile.write("\n")

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

def run_pling(file1,file2,iteration,parralel,i):
    """Runs pling 5 times on the fixed pair (file1, file2) — no fragment
    reordering — and appends the 5 resulting distances to
    allDistances_v2_base_pling.txt. `i` (the loop counter passed in) is
    shadowed by the internal `for i in range(5)` loop, so the caller's
    value of `i` is unused inside the function."""

    path = subprocess.check_output("echo $PWD",shell=True,text=True).strip()

    
    distances_sample = []

	#creating the input file for pling

    g = open(f"./py_out/input{iteration}_{parralel}.txt","w")
    g.write(path+file2+'\n')
    g.write(path+file1+"\n")
    g.close()

    for i in range(5):
        subprocess.run(
	    f"pling $PWD/py_out/input{iteration}_{parralel}.txt output_dir_try{i} align --containment_distance 1 --dcj 10",
		    		shell=True, check=True
		)
        distances_sample.append(int(distance(f"output_dir_try{i}")))
        subprocess.call(f"rm -r ./output_dir_try{i}", shell=True)
	
    print(len(distances_sample))
    w = open("allDistances_v2_base_pling.txt", "a")
    writeTo(w, distances_sample, file1, file2)
    w.close()
	
