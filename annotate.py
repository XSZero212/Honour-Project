"""
annotate.py — exploratory/dead code, kept for reference only.

The `comparisons`/`distances` module-level setup this file's plotting
functions depend on is commented out (lines 13-29 below), so none of
plot_distances()/plot_means()/plot_variances()/plot_mash_vs_pling() can run
as-is — they reference `distances`, `comparisons`, and `colors` globals
that are never defined in this file. There is also no `__main__` block, so
nothing here executes on its own; this was a scratch file for looking up
each pair's mash distance (via `egrep` against distances.tsv) alongside
its pling distance, before that logic moved into plot_vs_ground_trhuth.py.
Kept for reference rather than deleted, per project owner's request.
"""
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import sys
import pandas as pd

def writeTo(resultFile,distances):
    """Writes one value per line to resultFile."""
    g = open(resultFile,"w")
    for i in distances:
        g.write(str(i)+"\n")
    g.close()

# comparisons = []
# distances = []
# with open("all_pling_Distances.txt") as f:
# 	for line in f:
# 		comparisons.append((line.split(" ")[0].strip(),line.split(" ")[1].strip()))
# 		distances.append(line.split(" ")[2].strip())

# comparisons = comparisons[1:]
# distances = distances[1:]
# distances = np.array(distances)


# similarities = []
# for i in comparisons:
#  	similarities.append(subprocess.check_output(f'egrep "{i[0]}.*{i[1]}" ./distances.tsv| cut -d "\t" -f3 ',shell=True,text=True).strip())

# writeTo("similarities.txt",similarities)


def plot_distances():
	# Requires module-level `distances` (2D: one row per plasmid, values =
	# sampled distances) and `colors` — neither is defined in this file (see
	# module docstring). Plots each row's values against its own first
	# column (used as a stand-in "ground truth") to sanity-check spread.
	groundTruth =[]
	for i in range(len(distances)):
		groundTruth.append([distances[i,0] for x in distances[i]])
		plt.scatter(groundTruth[i],distances[i],color=colors[i][:len(distances)])
	plt.xlabel("Pling Distance")
	plt.ylabel("Permutated Distance")

	groundTruth = np.array(groundTruth)
	print(groundTruth.min())
	linex = [groundTruth.min(),groundTruth.max()]
	plt.plot(linex,linex)

def plot_means():
	# Same missing-globals caveat as plot_distances(). Plots each pair's
	# mean sampled distance against its first-column value, and prints the
	# MSE between them.
	groundTruth =[]
	means = []
	for i in range(len(distances)):
		mean = 0
		for j in distances[i,:]:
			mean+=j
		mean/=len(distances[i,:])
		groundTruth.append([distances[i,0] for x in distances[i]])
		means.append([mean for x in distances[i]])
		plt.scatter(groundTruth[i],means[i],color=colors[i][:len(distances)])
	plt.xlabel("Pling Distance")
	plt.ylabel("Mean")
	groundTruth = np.array(groundTruth)
	means  = np.array(means)
	print(groundTruth.min())
	linex = [groundTruth.min(),groundTruth.max()]

	print(np.square(np.subtract(groundTruth[:,0],means[:,0])).mean())
	plt.plot(linex,linex)

def plot_variances():
	# Requires module-level `distances` and `comparisons` (see module
	# docstring). For each pair, computes the variance of its sampled
	# distances, looks up the corresponding mash distance via `egrep`
	# against distances.tsv, and correlates/plots variance vs mash distance.
	variances = []
	for i in distances:
		mean = 0
		for j in i:
			mean += int(j)
		mean/=len(i)
		vari = 0
		for j in i:
			vari += (int(j)-mean)**2
		variances.append(float(vari/(len(i)-1)))
	similarities = []
	for i in comparisons:
		similarities.append(float(subprocess.check_output(f'egrep "{i[0]}.*{i[1]}" ./distances.tsv| tr "\t" " " | cut -d " " -f3 ',shell=True,text=True).strip()))

	print(np.corrcoef(variances,similarities))
	plt.scatter(similarities,variances)
	plt.xlabel("Mash Distance")
	plt.ylabel("Variances")

def plot_mash_vs_pling():
	# Requires module-level `distances` and `comparisons`. Plots each
	# pair's first sampled distance value against its mash distance
	# (looked up via `egrep` against distances.tsv).
	distance = [i[0] for i in distances]
	similarities = []
	for i in comparisons:
		similarities.append(float(subprocess.check_output(f'egrep "{i[0]}.*{i[1]}" ./distances.tsv| tr "\t" " " | cut -d " " -f3 ',shell=True,text=True).strip()))

	print(np.corrcoef(distance,similarities));
	plt.scatter(similarities,distance)
	plt.xlabel("Mash Distance")
	plt.ylabel("Pling Distance")
