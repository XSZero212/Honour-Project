"""
break_distance.py — this project's own alignment-based rearrangement
distance, as an alternative to calling out to pling.

Pipeline (entry point: my_implementation(plasmidA, plasmidB)):
  1. Align plasmidA against plasmidB with nucmer/delta-filter/show-coords.
  2. cluster_alignments(): greedily pick the best non-overlapping alignment
     blocks on each side (by length, then identity), then run CD-HIT on the
     aligned regions of each plasmid to cluster near-identical blocks into
     shared "synteny block" IDs — this is what turns two plasmid sequences
     into the block-permutation representation (a "genome" as used by the
     2-break-distance code below: a list of signed block IDs, sign
     encoding orientation) that the double-cut-and-join / 2-break distance
     formula operates on.
  3. calculateGaps(): penalises regions of each plasmid that weren't
     covered by any selected alignment block (unaligned "indels"), added on
     top of the 2-break distance.
  4. calculate2BreakDistance(): classic genome-rearrangement 2-break
     distance via the "colored edges" breakpoint-graph construction —
     builds the graph on both genomes' edges, union-finds connected
     components (cycles), and returns (#blocks - #cycles).

This is a from-scratch, alignment-driven distance measure, distinct from
pling's own algorithm; generateShortRead.py can be pointed at either
implementation via its `implementation` flag.
"""
import re
import subprocess
import math
import matplotlib.pyplot as plt
import copy

def overlaps(aln,used):
    """True if the alignment interval `aln` (a show-coords line for one
    side of an alignment) is fully contained within any interval already
    in `used` — used to skip redundant/nested alignments when greedily
    selecting non-overlapping blocks in cluster_alignments()."""
    distance_p = re.split(" +",aln)
    for alg in used:
        distance_2 = re.split(" +",alg)
        if int(distance_p[1])>=int(distance_2[1]) and int(distance_p[2])<=int(distance_2[2]):
            return True
    return False

def writeForClustering(plasmid,used):
    """Extracts the sequence under each selected alignment interval in
    `used` from `plasmid`'s FASTA, writes them as individually-numbered
    records to cluster_plasmid.fna, then runs CD-HIT-EST (95% identity) to
    group near-identical blocks together. Returns (assignedClusters, blocks)
    where assignedClusters maps a block's numeric ID to its CD-HIT cluster
    ID, and blocks is the list of (id, start, end) intervals — the cluster
    ID is what becomes the shared synteny-block label between plasmidA and
    plasmidB in cluster_alignments()."""
    #writes synteny blocks to a file so we can run CD-HIT on them
    f = open("cluster_plasmid.fna","w+")
    g = open(f"./{plasmid}","r")
    g.readline()
    text = g.read().replace("\n","")
    g.close()

    newPlasmid = ""
    i = 1
    blocks = []
    for block in used:
        coords = re.split(" +",block)
        
        if int(coords[1])<int(coords[2]):
            newPlasmid+=f">{i}\n"
            newPlasmid+=text[int(coords[1]):int(coords[2])+1]
            newPlasmid+="\n"
            blocks.append((i,int(coords[1]),int(coords[2])))
        else:
            newPlasmid+=f">{i}\n"
            newPlasmid+=text[int(coords[2]):int(coords[1])+1]
            blocks.append((i,int(coords[1]),int(coords[2])))
            newPlasmid+="\n"
        i+=1
    f.write(newPlasmid)
    f.close()

    #we have to run CD-HIT on them
    subprocess.call("cd-hit-est -i cluster_plasmid.fna -o clustered_blocks.fasta -c 0.95 -n 10",shell=True)

    currentCluster = 0
    assignedClusters = {}
    with open("clustered_blocks.fasta.clstr") as file:
        for line in file:
            if "Cluster" in line:
                currentCluster+=1
            else:
                plasmidId = int(line.split(">")[1].split("...")[0])
                assignedClusters[plasmidId] = currentCluster

    return assignedClusters,blocks

def find_cluster_id(assignedClusters,blocks,distance):
    """Looks up which CD-HIT cluster the alignment interval `distance`
    (a parsed show-coords start/end pair) belongs to, by matching its
    coordinates against the (id, start, end) tuples in `blocks`."""
    for block in blocks:
        if block[1]==int(distance[1]) and block[2]==int(distance[2]):
            return assignedClusters[block[0]]

def removeDuplicates(blocks_p,blocks_q):
    """Drops repeated cluster-ID occurrences within each block sequence
    (keeping the first occurrence, by absolute value/orientation-agnostic),
    since the 2-break distance formula expects each synteny block to
    appear exactly once per genome. `copies` counts how many were dropped
    (used as part of the gap/indel penalty in run())."""
    used_p = set()
    used_q = set()
    result_p = []
    result_q = []
    copies = 0
    for id,pos in enumerate(blocks_p):
        if abs(pos) not in used_p:
            used_p.add(abs(pos))
            result_p.append(blocks_p[id])
        else:
            copies+=1
        
    for id,pos in enumerate(blocks_q):
        if abs(pos) not in used_q:
            used_q.add(abs(pos))
            result_q.append(blocks_q[id])
        copies+=1
   
    return [result_p],[result_q],copies

def calculateGaps(blocks,plasmid):
    """Penalises unaligned regions of `plasmid` (gaps of >= 500bp between
    consecutive selected blocks, or before the first / after the last
    block). Returns (1 - total_unaligned_bp/plasmid_length) * num_gaps —
    a fraction-of-plasmid-aligned term scaled by how many separate gaps
    there are, used by run() as part of the total distance."""
    g = open(f"./{plasmid}","r")
    g.readline()
    text = g.read().replace("\n","")
    g.close()
    notOg = copy.deepcopy(blocks)
    notOg = sorted(notOg,key=lambda x:x[0])
    thresh = 500
    gaps = 0
    totalDist = 0
    if notOg[0][0]>thresh:
        gaps+=1
        totalDist+=notOg[0][0]
    for i in range(len(notOg)-1):
        if notOg[i+1][0]-notOg[i][1]>=thresh:
            gaps+=1
            totalDist+=notOg[i+1][0]-notOg[i][1]
    if blocks[-1][1]<len(text)-thresh:
        gaps+=1
        totalDist+=len(text)-blocks[-1][1]

    return (1-totalDist/len(text))*gaps

def cluster_alignments(plasmidA,plasmidB,file1):
    """Builds the shared synteny-block representation of plasmidA and
    plasmidB from a show-coords alignment file (`file1`).

    Parses each alignment line into (intervalA, intervalB, length,
    identity), sorts by length then identity (best first), and greedily
    keeps only alignments whose A- and B-side intervals don't overlap any
    already-selected block (overlaps()). The kept intervals from each side
    are separately clustered via CD-HIT (writeForClustering()) to assign
    matching block IDs across the two plasmids, then encoded as signed
    integer sequences (sign = strand orientation) — the same "genome"
    format used by calculate2BreakDistance().

    Returns ((genomeA, genomeB, copies), unaligned) where genomeA/genomeB
    are the block-ID sequences with duplicates removed (removeDuplicates()),
    copies is the number of duplicate blocks dropped, and unaligned is the
    combined gap penalty from calculateGaps() on both plasmids."""
    blocks_p = []
    blocks_q = []
    info = []
    with open(file1) as g:
        for line in g:
            if "|" in line and "[S1]" not in line:
                info.append((line.split("|")[0],line.split("|")[1],int(re.split(" +",line.split("|")[2])[1]),float(re.split(" +",line.split("|")[3])[1])))
                
    

    info = sorted(info,key=lambda x: (-x[2],-x[3]))

    selected = []
    usedA = []
    usedB = []


    for aln in info:
        if not overlaps(aln[0], usedA) and not overlaps(aln[1], usedB):
            selected.append(aln)
            usedA.append(aln[0])
            usedB.append(aln[1])
                
    assignedClustersA,blocksA = writeForClustering(plasmidA,usedA)
    assignedClustersB,blocksB = writeForClustering(plasmidB,usedB)

    for aln in selected:
       
        distance_p = re.split(" +",aln[0])
        distance_q = re.split(" +",aln[1])

        clusterId_p = find_cluster_id(assignedClustersA,blocksA,distance_p)
        clusterId_q = find_cluster_id(assignedClustersB,blocksB,distance_q)

        if int(distance_p[1])>int(distance_p[2]):
            blocks_p.append((int(distance_p[2]),int(distance_p[1]),-1*clusterId_p))
        else:
            blocks_p.append((int(distance_p[1]),int(distance_p[2]),clusterId_p))
                
        if int(distance_q[1])>int(distance_q[2]):
            blocks_q.append((int(distance_q[2]),int(distance_q[1]),-1*clusterId_q))
        else:
            blocks_q.append((int(distance_q[1]),int(distance_q[2]),clusterId_q))
                
    
    blocks_q = sorted(blocks_q, key=lambda tup: tup[1])
    unalignedP = calculateGaps(blocks_p,plasmidA)
    unalignedQ = calculateGaps(blocks_q,plasmidB)

    return removeDuplicates([tup[2] for tup in blocks_p],[tup[2] for tup in blocks_q]),unalignedP+unalignedQ

def read_alignemtns(file1):
    """Simpler alternative to cluster_alignments(): treats every alignment
    line in the show-coords file `file1` as its own synteny block (numbered
    in file order, no CD-HIT clustering, no overlap filtering). Not called
    by run()/my_implementation() — kept as a lighter-weight variant."""
    i=1
    blocks_p = []
    blocks_q = []
    with open(file1) as g:
        for line in g:
            if "|" in line and "[S1]" not in line:
                distances_big = line.split("|")[0:2]
                distance_p = re.split(" +",distances_big[0])
                distance_q = re.split(" +",distances_big[1])
                
                if int(distance_p[1])>int(distance_p[2]):
                    blocks_p.append((int(distance_p[2]),int(distance_p[1]),-i))
                else:
                    blocks_p.append((int(distance_p[1]),int(distance_p[2]),i))
                
                if int(distance_q[1])>int(distance_q[2]):
                    blocks_q.append((int(distance_q[2]),int(distance_q[1]),-i))
                else:
                    blocks_q.append((int(distance_q[1]),int(distance_q[2]),i))
                
                i+=1

    blocks_q = sorted(blocks_q, key=lambda tup: tup[1])

    return [[tup[2] for tup in blocks_p]],[[tup[2] for tup in blocks_q]]


def run(plasmidA,plasmidB):
    """Combines the synteny-block distance and the gap/indel penalty into
    one score: penal * (unaligned_penalty + duplicate_block_count) +
    2-break distance on the deduplicated block sequences. Assumes
    plasmid.coords (produced by my_implementation()) is the current
    nucmer/show-coords alignment between plasmidA and plasmidB."""
    genomes, indel_distance = cluster_alignments(plasmidA,plasmidB,"plasmid.coords")
    # print(genomes_p)
    # print()
    # print(genomes_q)
    # print()
    penal=1
    dist = penal*(indel_distance+genomes[2])+calculate2BreakDistance(genomes[0], genomes[1])
    return dist
    
def readGenomeFromFile(file):
        """Parses the "(±1 ±2 ...)(±3 ...)..." chromosome-string format
        (as written by syntethic.write_test()) into the list-of-signed-int
        genome representation used throughout this module. Used only by
        test()."""
        f = open(f'{file}', 'r')
        data = []
        for line in f:
            data.append(line.strip())

        for g in data:
            g = g.split(')(')
            genome = []
            for d in g:
                d = d[:-1].split()
                if d[0][0]=="(":
                    d[0]=d[0][1:]
                genome.append([int(i[1:]) if i[0]=="+" else -1*int(i[1:]) for i in d])
                # genome.append([int(d[0][1:] if '('==d[0][0] else d[0])] + [int(e) for e in d[1:-1]] +\
                # [int(d[-1][:-1] if ')'==d[-1][-1] else d[-1])])
            return genome


def chromosomeToCycle(chromosome):
    """Standard genome-rearrangement construction: expands each signed
    block in `chromosome` into its two graph nodes (head/tail), oriented
    according to the block's sign, per Bafna-Pevzner / Compeau-Pevzner
    2-break-distance theory."""
    l = len(chromosome)
    nodes = [0]*(2*l)
    for j in range(l):
        i = chromosome[j]
        if i > 0:
            nodes[2*j] = 2*i-1
            nodes[2*j+1] = 2*i
        else:
            nodes[2*j] = -2*i
            nodes[2*j+1] = -2*i-1
    return nodes
        
def cycleToChromosome(nodes):
    """Inverse of chromosomeToCycle(): converts a node-cycle back into a
    signed-block chromosome. Not called elsewhere in this module — kept
    for completeness/round-tripping."""
    l = len(nodes) // 2
    chromosome = [0]*l
    for j in range(l):
        if nodes[2*j] < nodes[2*j+1]:
            chromosome[j] = nodes[2*j+1]//2
        else:
            chromosome[j] = -nodes[2*j]//2
    return chromosome
    

def coloredEdges(genome):
    """Builds the set of "colored edges" (breakpoint-graph edges) for one
    genome — one color per genome, so calculate2BreakDistance() unions
    plasmidA's edges with plasmidB's edges to build the joint breakpoint
    graph whose cycle count gives the 2-break distance."""
    # transforms the geome that you give in a set of edges
    edges = set()
    for chromosome in genome:
        nodes = chromosomeToCycle(chromosome)
        nodes.append(nodes[0])
        for j in range(len(chromosome)):
            edges.add(((nodes[2*j+1], nodes[2*j+2])))
    return edges
        
def calculate2BreakDistance(P, Q):
    """Classic 2-break distance between two genomes P and Q (each a list of
    chromosomes, each chromosome a list of signed synteny-block IDs).
    Builds the joint breakpoint graph (coloredEdges(P) ∪ coloredEdges(Q)),
    finds its connected components via union-find, and returns
    (total_blocks - number_of_cycles) — the minimum number of 2-break
    (double-cut-and-join) operations needed to transform P into Q.
    (dist_formula_from_the_file / dist_updated are alternate odd-cycle-based
    formulas kept for reference/comparison; the function returns
    dist_2break_original.)"""
    blocks = sum([len(a) for a in P])
    edges = coloredEdges(P).union(coloredEdges(Q))
    parent = dict()
    rank = dict()
    for e in edges:
        parent[e[0]] = e[0]
        parent[e[1]] = e[1]
        rank[e[0]] = 0
        rank[e[1]] = 0

    def findParent(i):
        if i != parent[i]:
            parent[i] = findParent(parent[i])
        return parent[i]
        
    def union(i, j):
        i_id = findParent(i)
        j_id = findParent(j)
        if i_id == j_id:
            return
        if rank[i_id] > rank[j_id]:
            parent[j_id] = i_id
        else:
            parent[i_id] = j_id
            if rank[i_id] == rank[j_id]:
                rank[j_id] += 1
        
    for e in edges:
        union(e[0], e[1])

    cycles = dict()

    for e in edges:
        id = findParent(e[0])

        if id in cycles:
            cycles[id]+=1
        else:
            cycles[id]=1
        
    odd_cycles = 0
    for i in cycles:
        if (cycles[i]/2)%2==1:
            odd_cycles+=1
    dist_formula_from_the_file = math.ceil((blocks - odd_cycles)/2)
    dist_updated = math.ceil((blocks - len(cycles))/2)
    dist_2break_original = (blocks - len(cycles))
    return dist_2break_original

def printGenome(genome):
        result = ''
        for chromosome in genome:
            result += '('+' '.join(['+'+str(e) if e>0 else str(e) for e in chromosome])+')'
        print(result)

def my_implementation(plasmidA,plasmidB):
    """Public entry point used by generateShortRead.py when
    implementation=True. Aligns plasmidA against plasmidB with nucmer
    (writing out.delta -> out.filtered.delta -> plasmid.coords), then
    calls run() to turn that alignment into a rearrangement distance."""
    subprocess.call(f"nucmer --mum --prefix=out --diagdiff=20 --breaklen=50 --minmatch 50 ./{plasmidA} ./{plasmidB} ",shell=True)
    subprocess.call("delta-filter -1 out.delta > out.filtered.delta",shell=True)
    subprocess.call("show-coords -rcl out.filtered.delta > plasmid.coords",shell=True)
    return run(plasmidA,plasmidB)

def test():
    """Sanity check for calculate2BreakDistance() against synthetic data
    with a known ground-truth number of rearrangement events (encoded in
    the filenames under syntethic/fastas, as written by syntethic.write_test()).
    Compares the computed 2-break distance to the known event count and
    scatter-plots ground truth vs computed distance."""
    plasmid_files = subprocess.check_output("ls -1 syntethic/fastas",shell=True,text=True)
    files = plasmid_files.split("\n")[:-1]
    same = 0
    total = 0
    not_same = []
   
    for f in files:
        text = f.split("fna")
        plasmid = text[0]+"fna"
        changes = int(text[1].split("_")[1])
        nrOfFragments = int(text[1].split("_")[2])
        genomep = [[i for i in range(nrOfFragments)]]
        genomeq = readGenomeFromFile(f"syntethic/fastas/{f}")
        dist = calculate2BreakDistance(genomep, genomeq)
        strange = []
        not_same.append((changes,dist))
        if dist == changes:
            same+=1
        elif dist>changes:
            strange.append((changes,dist))
        
        total+=1
    print(same)
    print(same/total)
    print(not_same)
    print(strange)
    plt.scatter([tup[0] for tup in not_same],[tup[0] for tup in not_same],color="red")
    plt.scatter([tup[0] for tup in not_same],[tup[1] for tup in not_same])
    plt.xlabel("The ground truth distance")
    plt.ylabel("The generated distance")
    plt.title("2-break distance")
    plt.show()

#   nr_of_changes = []
#     mean = 0
#     for f in files:
#         changes = int(f.split("fna")[1].split("_")[1])
#         nr_of_changes.append(changes)
#         mean+=changes
#     print(mean/len(nr_of_changes))
#     # plt.plot(nr_of_changes)
#     # plt.show()

def graph_construction(genomeP,genomeQ):
    """Alternative single-chromosome 2-break-distance implementation: builds
    the breakpoint graph directly as black edges (genomeP's adjacencies)
    and grey edges (genomeQ's adjacencies), union-finds cycles while
    tracking black-edge parity per component, and returns
    (num_blocks - odd_cycles) / 2. Not called elsewhere in this module —
    kept as an alternate/earlier formula next to calculate2BreakDistance()."""
    nodes = [0]*(2*len(genomeP))
    edges = []

    #generating the black edges
    for i in range(len(genomeP[0])):
        idx = (i+1)%len(genomeP[0])
        head = 2*i+1 if genomeP[0][i]>=0 else 2*i
        tail = 2*idx if genomeP[0][i]>=0 else 2*idx+1
        edges.append(((head,tail),0))
    #generating the grey edges
    for i in range(len(genomeQ[0])):
        idx = (i+1)%len(genomeQ[0])
        head = 2*i+1 if genomeQ[0][i]>=0 else 2*i
        tail = 2*idx if genomeQ[0][i]>=0 else 2*idx+1
        edges.append(((head,tail),1))
    
    cycles = dict()
    nrOfCycles = 0

    parent = dict()
    rank = dict()
    black = dict()
    for e in nodes:
        parent[e] = e
        rank[e] = 0
        black[e] = 0

    def findParent(i):
        if i != parent[i]:
            parent[i] = findParent(parent[i])
        return parent[i]
        
 
    #the actual calculation
    for i in edges:
        parent_head = findParent(i[0][0])
        parent_tail = findParent(i[0][1])
        if parent_tail!=parent_head:
            black[parent_head]+=black[parent_tail]
            parent[parent_tail]=parent_head
            if i[1]==0:
                black[parent_head]+=1

    nodesSets = set()

    for e in edges:
        id = findParent(e[0][0])
        nodesSets.add(id)
    
    oddcycles = 0
    for node in nodesSets:
        if black[node]%2==1:
            oddcycles+=1
    
    return math.ceil(len(genomeP[0])-oddcycles)/2
        
#test()
            


