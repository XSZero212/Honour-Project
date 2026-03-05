import subprocess
import matplotlib.pyplot as plt
import numpy as np
import random as rm
import sys
from syntethic import run
#splits the distances into four groups


def create_input():
	syntethic_fasta = subprocess.check_output("ls -1 ./syntethic/regenerated/",shell=True,text=True)
	syntethic_fasta = syntethic_fasta.split("\n")
	f = open("input_syntethic.txt","w+")
	for file in syntethic_fasta:
		og_fasta = file.split("_")[0]
		f.write(f"./fastas/{og_fasta} ./syntethic/regenerated/{file}\n")
	f.close()

i=0
def sample_random_plasmids():
	pair_plasmids = []
	with open("uniqueDistancesv2.txt") as f:
		for x in f:
			s = x.strip()
			pair_plasmids.append((s.split()[0][1:],s.split()[1][1:]))

	i=0
	while True:
		plasmids_pairs = rm.sample(pair_plasmids,10)
		commands = [f"python ./generateShortRead.py {j[0]} {j[1]} {i} {index} {4}" for index,j in enumerate(plasmids_pairs[:])]
		procs = [subprocess.Popen(j,shell=True) for j in commands]
		for p in procs:
			p.wait()
		print(f"The batch finished {i}")
		i+=1

def sample_continous(implementation):
	pair_plasmids = []
	with open("input_syntethic.txt") as f:
		for x in f:
			s = x.strip()
			pair_plasmids.append((s.split()[0][1:],s.split()[1][1:]))

	i=0
	k=0
	parralel = 10
	while i<len(pair_plasmids):
		commands = [f"python ./generateShortRead.py {j[0]} {j[1]} {k} {index} {6} {implementation}" for index,j in enumerate(pair_plasmids[i:i+parralel])]
		procs = [subprocess.Popen(j,shell=True) for j in commands]
		for p in procs:
			p.wait()
		if not implementation:
			for index,_ in enumerate(pair_plasmids[i:i+parralel]):
				subprocess.call(f"rm -r ./py_out/input{k}_{index}.txt", shell=True)
				subprocess.call(f"rm -r ./py_out/try{k}_{index}.fasta", shell=True)
		print(f"The batch finished {i}")
		i+=parralel
		k+=1

def change_fragments_known(size,implementation):
	plasmid_pairs = []
	with open("allDistances_v2_1.txt") as f:
		for line in f:
			if "/" in line:
				plasmid_pairs.append(("/fastas/"+line.split("/")[2].strip(),"/fastas/"+line.split("/")[4].strip()))
	i=0
	k=0
	parralel = 1
	while i<len(plasmid_pairs):
		commands = [f"python ./generateShortRead.py {j[0]} {j[1]} {k} {index} {size} {implementation}" for index,j in enumerate(plasmid_pairs[i:i+parralel])]
		procs = [subprocess.Popen(j,shell=True) for j in commands]
		for p in procs:
			p.wait()
		if not implementation:
			for index,_ in enumerate(plasmid_pairs[i:i+parralel]):
				subprocess.call(f"rm -r ./py_out/input{k}_{index}.txt", shell=True)
				subprocess.call(f"rm -r ./py_out/try{k}_{index}.fasta", shell=True)
		print(f"The batch finished {i}")
		i+=parralel
		k+=1
#print(commands)
		#subprocess.call(f"python ./generateShortRead.py {s.split()[0][1:]} {s.split()[1][1:]} {i}", shell=True)

def trying_something():
	plasmid_pairs = []
	with open("allDistances_v2_5.txt") as f:
		for line in f:
			if "/" in line:
				plasmid_pairs.append(("/fastas/"+line.split("/")[2].strip(),"/fastas/"+line.split("/")[4].strip()))
	
	subprocess.call(f"python ./generateShortRead.py {plasmid_pairs[0][0]} {plasmid_pairs[0][1]} {1} {0} {5}",shell=True)

def generate_syntethic():
	plasmid_pairs = set()
	with open("allDistances_v2_1.txt") as f:
		for line in f:
			if "/" in line:
				plasmid_pairs.add("/fastas/"+line.split("/")[2].strip())
				plasmid_pairs.add("/fastas/"+line.split("/")[4].strip())
	k=30
	plasmid = rm.sample(sorted(plasmid_pairs),k)
	for i in range(k):
		nr_of_fragments = rm.randrange(1,20)
		run(plasmid[i],nr_of_fragments)

def run_syntethic(size):
	plasmid_pairs = set()
	plasmid_files = subprocess.check_output("ls -1 syntethic/regenerated",shell=True,text=True)
	files = plasmid_files.split("\n")[:-1]
	
	i=0
	for f in files:
		og = f.split(".fna")[0]
		subprocess.call(f"python ./generateShortRead.py /fastas/{og}.fna /syntethic/regenerated/{f} {i} {1} {size} True",shell=True)
		i+=1

if __name__=="__main__":
	func = int(sys.argv[1])
	size = int(sys.argv[2])
	implementation = sys.argv[3].lower()== "true"
	if func == 0:
		sample_random_plasmids()
	elif func==1:
		change_fragments_known(size,implementation)
	elif func==2:
		generate_syntethic()
	elif func==3:
		run_syntethic(size)
	elif func==4:
		sample_continous(implementation)
	else:
		trying_something()
