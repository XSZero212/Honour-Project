import subprocess
import matplotlib.pyplot as plt
import numpy as np
import sys
import pandas as pd
import scipy.stats as stats


def short_write():
	distances =[]
	comparisons = []
	with open("allDistances_v2_3.txt") as f:
		for line in f:
			if "/" in line:
				comparisons.append((line.split("/")[2].replace(".1.fna",""),line.split("/")[4].strip().replace(".1.fna","")))
			else:
				distances.append(np.fromstring(line.strip(),sep=" "))
	h = open("ground_Thruts.txt","w")
	for i in range(len(comparisons)):
		h.write(comparisons[i][0]+" "+comparisons[i][1]+" "+str(distances[i][0])+"\n")
	h.close()

def read(file):
	distances = []
	comparisons = []
	groundTruth = []
	with open(file) as f:
		for line in f:
			if "/" in line:
				comparisons.append((line.split("/")[2].replace(".1.fna",""),line.split("/")[4].strip().replace(".1.fna","")))
				groundTruth.append(float(subprocess.check_output(f'egrep "{line.split("/")[2].replace(".1.fna","")} {line.split("/")[4].strip().replace(".1.fna","")}" ./ground_Thruts.txt|cut -d " " -f3 ',shell=True,text=True).strip()))
			else:
				distances.append(np.fromstring(line.strip(),sep=" "))
	
	return np.array(comparisons),np.array(distances),np.array(groundTruth)
def short_read(file):
	distances = []
	comparisons = []
	groundTruth = []
	with open(file) as f:
		for line in f:
			if "/" in line:
				comparisons.append((line.split("/")[2].replace(".1.fna",""),line.split("/")[5].strip()))
				fna = line.split("/")[5].strip().split(".fna")[1]
				groundTruth.append(float(fna.split("_")[1]))

			else:
				distances.append(np.fromstring(line.strip(),sep=" "))
	
	return np.array(comparisons),np.array(distances),np.array(groundTruth)

def plot_distances(distances,ground_Truth,axs,x,y,size):
	groundTruth =[]
	for i in range(len(distances)):
		groundTruth.append([ground_Truth[i] for _ in distances[i]])
		if x==0 and y==0:
			axs.scatter(groundTruth[i],distances[i])
		else:
			axs[x,y].scatter(groundTruth[i],distances[i])
	groundTruth = np.array(groundTruth)
	print(groundTruth.min())
	linex = [groundTruth.min(),groundTruth.max()]

	if x==0 and y==0:
		axs.set_title(f"Number Of Fragments: {size}")
		axs.plot(linex,linex)
		axs.set_xlabel("Ground Truth")
		axs.set_ylabel("Calculated Distance")
	else:
		axs[x,y].set_title(f"Number Of Fragments: {size}")
		axs[x,y].plot(linex,linex)

def plot_means(distances,ground_Truth,axs,x,y,size):
	groundTruth =[]
	means = []
	for i in range(len(distances)):
		mean = 0
		for j in distances[i,:]:
			mean+=j
		mean/=len(distances[i,:])
		groundTruth.append([ground_Truth[i] for _ in distances[i]])
		means.append([mean for x in distances[i]])
		if x==0 and y==0:
			axs.scatter(groundTruth[i],means[i])
		else:
			axs[x,y].scatter(groundTruth[i],means[i])

	groundTruth = np.array(groundTruth)
	means  = np.array(means)
	print(groundTruth.min())
	linex = [groundTruth.min(),groundTruth.max()]

	print(np.square(np.subtract(groundTruth[:,0],means[:,0])).mean())

	if x==0 and y==0:
		axs.set_title(f"Number Of Fragments: {size}")
		#axs.scatter(ground_Truth,ground_Truth,color="red")
		axs.plot(linex,linex)
		axs.set_xlabel("Ground Truth")
		axs.set_ylabel("Mean Distance")
	else:
		axs[x,y].set_title(f"Number Of Fragments: {size}")
		axs[x,y].plot(linex,linex)

def plot_means_v2(distances,ground_Truth,axs=None,x=0,y=0,size=5):
	groundTruth =[]
	means = []
	for i in range(len(distances)):
		mean = 0
		for j in distances[i,:]:
			mean+=j
		mean/=len(distances[i,:])
		means.append(mean)
	if axs==None:
		plt.title(f"Number Of Fragments: {size}")
		groundTruth = np.array(groundTruth)
		means  = np.array(means)
		means = (means-ground_Truth)**2
		print(np.corrcoef(means,ground_Truth))
		tuples = [(means[i],ground_Truth[i]) for i in range(len(ground_Truth))]
		tuples.sort(key=lambda tup: tup[1])
		residuals = [t[0] for t in tuples]
		groundTruth = [t[1] for t in tuples]
		plt.scatter(groundTruth,residuals)
		plt.xlabel("Ground Truth Distance")
		plt.ylabel("Residuals")


	# axs[x,y].set_title(f"Number Of Fragments: {size}")
	# groundTruth = np.array(groundTruth)
	# means  = np.array(means)
	# means = (means-ground_Truth)**2
	# print(np.corrcoef(means,ground_Truth))
	# tuples = [(means[i],ground_Truth[i]) for i in range(len(ground_Truth))]
	# tuples.sort(key=lambda tup: tup[1])
	# residuals = [t[0] for t in tuples]
	# groundTruth = [t[1] for t in tuples]
	# axs[x,y].scatter(groundTruth,residuals)
	# axs[x,y].set_xlabel("Ground Truth Distance")
	# axs[x,y].set_ylabel("Residuals")

	return np.sum(means),size

def calculate_means(distance):
	means = []
	for i in range(len(distance)):
		mean = 0
		for j in distance[i,:]:
			mean+=j
		mean/=len(distance[i,:])
		means.append(mean)
	return means
def test_means(distance_3,distance_4,distance_6,distance_10):

	means_3 = calculate_means(distance_3)
	means_4 = calculate_means(distance_4)
	means_6 = calculate_means(distance_6)
	means_10 = calculate_means(distance_10)

	# Calculate the group means
	group_means = []
	group_means.append(np.mean(means_3))
	group_means.append(np.mean(means_4))
	group_means.append(np.mean(means_6))
	group_means.append(np.mean(means_10))


	# Set the significance level
	alpha = 0.05


	# Perform one-way ANOVA using the F-test
	
	f_value, p_value = stats.friedmanchisquare(*[means_3[:-1],means_4[:-1],means_6[:-1],means_10])

	# Print the results
	print("F-value:", f_value)
	print("p-value:", p_value)

	# Set the significance level
	alpha = 0.05

	# Check the hypothesis
	if p_value < alpha:
		print("Reject the null hypothesis")
	else:
		print("Fail to reject the null hypothesis")


def plot_variances(distances,comparisons,axs,x,y,size):
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
	binwidth=8
	axs[x,y].hist(variances,bins="fd")
	axs[x,y].set_title(f"Number of fragments: {size}")
	return np.mean(variances),size

def sorted_means(distances,ground_Truth):
	means = calculate_means(distances)
	errors  = (means-ground_Truth)**2
	tuples = [(errors[i],ground_Truth[i]) for i in range(len(ground_Truth))]
	tuples.sort(key=lambda tup:tup[1])
	return tuples
def sorted_min(distances,ground_Truths):
	minVals = []
	for i in distances:
		minVals.append(min(i))
	errors = (minVals-ground_Truths)**2
	tuples = [(errors[i],ground_Truths[i]) for i in range(len(ground_Truths))]
	tuples.sort(key=lambda tup:tup[1])
	return tuples

def compare_residuals(tuples_mean,tuples_min):
	mean = np.sum(np.sum([i[0] for i in tuples_mean]))
	mini = np.sum(np.sum([i[0] for i in tuples_min]))
	print("Mean: ")
	print(mean)
	print("Mini: ")
	print(mini)
	return (mean,mini)


comparisons_normal,distances_normal,ground_truth = short_read("allDistances_v2_6_True_2_break.txt")
#comparisons_changed,distances_changed,ground_truth_changed = short_read("allDistances_v2_6_True.txt")
# fig,axs = plt.subplots()
# plot_distances(distances_changed,ground_truth_changed,axs,x=0,y=0,size=6)
# plt.show()

# fig,axs = plt.subplots()
# plot_means(distances_changed,ground_truth_changed,axs,x=0,y=0,size=6)
# plt.show()

# fig,axs = plt.subplots(1,2)
# print(plot_means_v2(distances_changed,ground_truth_changed,x=0,y=0,size=6))
print(plot_means_v2(distances_normal,ground_truth,x=0,y=0,size=6))
plt.show()

# tuples_means_3 = sorted_means(distances_3,groundTrunth_3)
# tuples_min_3 = sorted_min(distances_3,groundTrunth_3)
# tuples_means_4 = sorted_means(distances_4,groundTrunth_4)
# tuples_min_4 = sorted_min(distances_4,groundTrunth_4)
# tuples_means_5 = sorted_means(distances_5,groundTrunth_5)
# tuples_min_5 = sorted_min(distances_5,groundTrunth_5)
# tuples_means_6 = sorted_means(distances_6,groundTrunth_6)
# tuples_min_6 = sorted_min(distances_6,groundTrunth_6)
# tuples_means_10 = sorted_means(distances_10,groundTrunth_10)
# tuples_min_10 = sorted_min(distances_10,groundTrunth_10)



# axs[0,0].set_title("Number of fragments: 3")
# axs[0,0].plot([i[1] for i in tuples_means_3],[i[0] for i in tuples_means_3],color="blue")
# axs[0,0].plot([i[1] for i in tuples_min_3],[i[0] for i in tuples_min_3],color="orange")
# axs[0,1].set_title("Number of fragments: 4")
# axs[0,1].plot([i[1] for i in tuples_means_4],[i[0] for i in tuples_means_4],color="blue")
# axs[0,1].plot([i[1] for i in tuples_min_4],[i[0] for i in tuples_min_4],color="orange")
# axs[1,0].set_title("Number of fragments: 6")
# axs[1,0].plot([i[1] for i in tuples_means_6],[i[0] for i in tuples_means_6],color="blue")
# axs[1,0].plot([i[1] for i in tuples_min_6],[i[0] for i in tuples_min_6],color="orange")
# axs[1,1].set_title("Number of fragments: 10")
# axs[1,1].plot([i[1] for i in tuples_means_10],[i[0] for i in tuples_means_10],color="blue")
# axs[1,1].plot([i[1] for i in tuples_min_10],[i[0] for i in tuples_min_10],color="orange")

# labels = ["Residual using average distance", "Residual using minimum distance"]
# fig.supxlabel("Ground Truth")
# fig.supylabel("Residual")
# fig.legend(labels=labels,loc="upper right")
# plt.show()
# values = []
# values.append(compare_residuals(tuples_means_3,tuples_min_3))
# values.append(compare_residuals(tuples_means_4,tuples_min_4))
# values.append(compare_residuals(tuples_means_6,tuples_min_6))
# values.append(compare_residuals(tuples_means_10,tuples_min_10))
# values.append(compare_residuals(tuples_means_5,tuples_min_5))

# plt.plot([3,4,5,6,10],[tup[0] for tup in values],color="blue")
# plt.plot([3,4,5,6,10],[tup[1] for tup in values],color="orange")
# labels = ["Residual using average distance", "Residual using minimum distance"]
# plt.legend(labels=labels,loc="upper right")
# plt.xlabel("Fragments")
# plt.ylabel("Residuals")
# plt.show()
#test_means(distances_3,distances_4,distances_6,distances_10)
# elbow = []
# elbow.append(plot_means_v2(distances_3,groundTrunth_3,axs,0,0,size=3))
# elbow.append(plot_means_v2(distances_4,groundTrunth_4,axs,0,1,size=4))
# #elbow.append(plot_means_v2(distances_5,groundTrunth_5,axs,0,0,size=5))
# elbow.append(plot_means_v2(distances_6,groundTrunth_6,axs,1,0,size=6))
# elbow.append(plot_means_v2(distances_10,groundTrunth_10,axs,1,1,size=10))

# elbow = np.array(elbow)
# fig.supxlabel("Variance")
# fig.supylabel("Frequency")
# plt.show()
# f_value, p_value = stats.friedmanchisquare(*[elbow[0,:],elbow[1,:],elbow[2,:],elbow[3,:]])

# # Print the results
# print("F-value:", f_value)
# print("p-value:", p_value)

# # Set the significance level
# alpha = 0.05

# # Check the hypothesis
# if p_value < alpha:
# 	print("Reject the null hypothesis")
# else:
# 	print("Fail to reject the null hypothesis")
# fig.supxlabel("Ground Truth")
# fig.supylabel("Residuals")

# res = [t[0] for t in elbow]
# sizes = [t[1] for t in elbow]
# plt.plot(sizes,res)
# plt.scatter(sizes,res,s=100)
# plt.ylabel("Average Variance")
# plt.xlabel("Number of fragments")
# plt.show()
def test():
    plasmid_files = subprocess.check_output("ls -1 syntethic/fastas",shell=True,text=True)
    files = plasmid_files.split("\n")[:-1]
    nr_of_changes = []
    mean = 0
    for f in files:
        changes = int(f.split("fna")[1].split("_")[1])
        nr_of_changes.append(changes)
        mean+=changes
    print(mean/len(nr_of_changes))

