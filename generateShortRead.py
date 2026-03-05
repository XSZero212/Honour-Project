import subprocess
import matplotlib.pyplot as plt
import numpy as np
import sys
import random as rdm
from break_distance import my_implementation

def writeTo(resultFile,distances,file1,file2):
	resultFile.write(file1+file2+"\n")
	for i in distances:
		resultFile.write(str(i)+" ")
	resultFile.write("\n")

def generateShort(nr,file,nrOfFragments) :
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

def try_something(val):
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
	#good reminder I quess, not exactly necessary but in case pling does not calculate the distance we add -1
	f = open(s+"/all_plasmids_distances.tsv")
	f.readline()
	values = f.readline()
	if values != "":
		return values.split('\t')[2].strip()
	else:
		return -1

def construct(indexes,iteration,parralel):
	#This constructs the fasta file by splitting the plasmids into sections of 70 chars
	#val called here starts from 1, because in val[0] we have the fasta index of the plasmid
	s =""
	for i in indexes:
		s+=val[i+1]
	i = 0
	#Adding newlines to the string
	newS = ""
	while i < len(s):
		newS += s[i:i+70]+'\n'
		i+=70
	#writting it to the respective file
	g = open(f"./py_out/try{iteration}_{parralel}.fasta","w")
	g.write(val[0]+'\n'+newS+'\n')
	g.close()
	
def generate(k,size,iteration,used,elements,parralel):
	# this function generates the permutations and then calculates their pling distance
	if k==size:
		global distances
		global var
		var+=1
		construct(elements,iteration,parralel)
		subprocess.call(f"pling $PWD/py_out/input{iteration}_{parralel}.txt output_dir_{var}_{parralel}_{size} align pling ./similarities.txt output align --containment_distance 1 --dcj 10 ",shell=True)
		distances.append(int(distance(f"output_dir_{var}_{parralel}_{size}")))
		subprocess.call(f"rm -r ./output_dir_{var}_{parralel}_{size}", shell=True)
	else:
		for i in range(size):
			if used[i]==False:
				elements[k]=i;
				used[i]=True;
				generate(k+1,size,iteration,used,elements,parralel);
				used[i]=False;

def factorial(k):
	prod = 1
	for i in range(1,k+1):
		prod*=i
	return prod
def reverse(n,size):
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

def empty(used):
	for i in range(1,len(used)):
		used[i] = False
	return used  
def fill(elems, k,used):
	for i in range(1,len(used)):
		if used[i]==False:
				elems[k] = i
				k+=1
	return elems
def get(n,used):
	k = 0;
	for i in range(1,len(used)):
		if used[i]==False:
			k+=1
		if k == n + 1:
			return i
	return -1
def nextSmall(used):
	for i in range(1,len(used)):
		if used[i]==False:
				return i;
	return len(used) - 1;

def subsample(size,file1,file2, iteration,parralel,implementation=False):
	distances_sample = []
	fact = factorial(size)
	for i in range(120):
		index = rdm.randint(0,fact-1)
		elements = reverse(index,size)
		elements = [i-1 for i in elements]
		construct(elements,iteration,parralel)

		if implementation:
			distances_sample.append(my_implementation(file2,f"./py_out/try{iteration}_{parralel}.fasta"))
		else:
		# I put parralel so different threads do not fight between themselves
		# I put sizes so that same parralel counts for different fragment sized do not fight amongst themselves
		# I put i so that same parralel calls do not fight against themselves
			subprocess.call(f"pling $PWD/py_out/input{iteration}_{parralel}.txt output_dir_{i}_{parralel}_{size} align --containment_distance 1 --dcj 10",shell=True)
			distances_sample.append(int(distance(f"output_dir_{i}_{parralel}_{size}")))
			subprocess.call(f"rm -r ./output_dir_{i}_{parralel}_{size}", shell=True)

	w = open(f"allDistances_v2_{size}_{implementation}.txt","a")
	writeTo(w,distances_sample,file1,file2)
	w.close()
	

def initial(file1,file2,iteration,parralel,nrOfFragments,implementation):
	
	#The number of lines of the fasta
	nr = subprocess.check_output(f"wc -l .{file1}",shell=True,text=True)
	nr = nr.strip().split(" ")[0]

	path = subprocess.check_output("echo $PWD",shell=True,text=True).strip()

	#this generates a list of plasmid fragments
	global val
	val = generateShort(int(nr),"."+file1,nrOfFragments)

	global distances
	distances = []

	#creating the input file for pling
	size = len(val)-1 
	used = [False]*size
	elements = [0]*size
	g = open(f"./py_out/input{iteration}_{parralel}.txt","w")
	g.write(path+file2+'\n')
	g.write(path+f"/py_out/try{iteration}_{parralel}.fasta")
	g.close()

	#calling the permutation algorithm
	if size<=5:
		global var
		var = 0
		generate(0,len(val)-1,iteration,used,elements,parralel)
		#saving the distances
		w = open(f"allDistances_v2_{size}.txt","a")
		writeTo(w,distances,file1,file2)
		w.close() 
	else:
		subsample(nrOfFragments,file1,file2,iteration=iteration,parralel=parralel)
	
def initial_my_implementation(file1,file2,iteration,parralel,nrOfFragments):
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


	