"""
syntethic.py — generates the synthetic plasmid database (syntethic/regenerated/)
by applying random rearrangement events to a real plasmid from fastas/, so
that the rearrangement distance between original and rearranged is known
exactly (ground truth for validating pling / break_distance.py against).

Entry point: run(file, nrOfFragments). Splits the original plasmid into
nrOfFragments equal blocks (split()), representing it as a single
chromosome [0, 1, ..., nrOfFragments-1], then for k=10 independent runs
applies a random number of REVERSAL/TRANSPOSITION/INSERTION/DELETION
events (generate()), each time re-measuring the true 2-break distance from
the original with TwoBreakDistance (so the reported ground truth accounts
for events that cancel out or produce an already-seen genome, rather than
just counting attempted events). Writes both a block-string representation
(write_test(), consumed by break_distance.readGenomeFromFile()) and the
reconstructed FASTA (write_fasta()) to syntethic/regenerated/, named
{plasmid}_{true_2break_distance}_{nrOfFragments}.fna.

TwoBreakDistance duplicates the 2-break-distance algorithm also implemented
in break_distance.py (calculate2BreakDistance/coloredEdges); kept here as a
self-contained class so this module doesn't depend on break_distance.py.
"""
from enum import Enum
import random
import subprocess
import copy
import matplotlib.pyplot as plt
import numpy as np

class Events(Enum):
    REVERSAL = 0
    TRANSLOCATION = 1
    TRANSPOSITION = 2
    INSERTION = 3
    DELETION = 4

class TwoBreakDistance:

    def __init__(self,plasmid1,plasmid2):
        #genomes = self.readGenomesFromFile()
        self.dist = self.calculate2BreakDistance(plasmid1, plasmid2)

    def chromosomeToCycle(self, chromosome):
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
        
    def cycleToChromosome(self, nodes):
        l = len(nodes) // 2
        chromosome = [0]*l
        for j in range(l):
            if nodes[2*j] < nodes[2*j+1]:
                chromosome[j] = nodes[2*j+1]//2
            else:
                chromosome[j] = -nodes[2*j]//2
        return chromosome
    
    def printChromosome(self, chromosome):
        print('('+' '.join(['+'+str(e) if e>0 else str(e) for e in chromosome])+')')

    def coloredEdges(self, genome):
        edges = set()
        for chromosome in genome:
            nodes = self.chromosomeToCycle(chromosome)
            nodes.append(nodes[0])
            for j in range(len(chromosome)):
                edges.add((nodes[2*j+1], nodes[2*j+2]))
        return edges
        
    def calculate2BreakDistance(self, P, Q):
        blocks = sum([len(a) for a in P])
        edges = self.coloredEdges(P).union(self.coloredEdges(Q))
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

        nodesSets = set()

        for e in edges:
            id = findParent(e[0])
            nodesSets.add(id)
        
        cycles = len(nodesSets)
        dist = blocks - cycles
        return dist

    def printGenome(self, genome):
        result = ''
        for chromosome in genome:
            result += '('+' '.join(['+'+str(e) if e>0 else str(e) for e in chromosome])+')'
        print(result)

def reversal(genome):
    """Reverses (and flips the sign/orientation of) a random contiguous
    segment within a random chromosome of `genome`, in place."""
    chr = random.randrange(len(genome))

    if len(genome[chr])<2:
        return genome

    i= random.randrange(len(genome[chr])-1)
    j = random.randrange(i+1,len(genome[chr]))

    segment = genome[chr][i:j+1]
    segment = [-x for x in reversed(segment)]

    genome[chr][i:j+1] = segment
    
    return genome

def translocation(genome):
    """Swaps a random suffix between two randomly chosen chromosomes.
                    Not used"""
    if len(genome) < 2:
        return genome

    i, j = random.sample(range(len(genome)), 2)
    chr1 = genome[i]
    chr2 = genome[j]

    if len(chr1) < 2 or len(chr2) < 2:
        return genome

    cut1 = random.randrange(1, len(chr1))
    cut2 = random.randrange(1, len(chr2))

    new_chr1 = chr1[:cut1] + chr2[cut2:]
    new_chr2 = chr2[:cut2] + chr1[cut1:]

    new_genome = genome.copy()
    new_genome[i] = new_chr1
    new_genome[j] = new_chr2

    return new_genome

def transposition(genome):
    """Swaps two whole chromosomes' positions in `genome`. As with
    translocation(), only meaningful for multi-chromosome genomes; plasmids
    here are modelled as a single chromosome."""
    if len(genome)<2:
        return genome
    
    i,j = random.sample(range(len(genome)),2)

    genome_inv = genome[i]
    genome[i] = genome[j]
    genome[j] = genome_inv

    return genome

def insertion(genome):
    """Duplicates a random block within chromosome 0 and inserts the copy
    at a random position (a "gene duplication" event, not sequence
    insertion from outside the genome)."""
    gene = random.randrange(len(genome[0]))
    position = random.randrange(len(genome[0]))

    genome[0].insert(position,genome[0][gene])

    return genome

def deletion(genome):
    """Removes a random block from chromosome 0 (no-op if only one block
    remains, since an empty chromosome isn't a valid genome)."""
    if len(genome[0])<=1 :
        return genome
    
    pos = random.randrange(len(genome[0]))
    genome[0].pop(pos)

    return genome

def reverse_complement(seq):
    complement = str.maketrans("ACGT", "TGCA")
    return seq.translate(complement)[::-1]

def write_fasta(plasmid_1,block_library,file,changes,nrOfFragments):
    """Reconstructs a FASTA sequence from the rearranged block-ID genome
    `plasmid_1` by concatenating each block's sequence from `block_library`
    (reverse-complementing blocks with a negative sign), and writes it to
    syntethic/regenerated/{plasmid}_{changes}_{nrOfFragments}.fna. `changes`
    is the true 2-break distance from the original, used to name the file
    with its own ground truth."""
    sequence = ""
    text = file.split("/")
    plasmid = text[2]
    sequence+=block_library[0]+'\n'
    for genome in plasmid_1:
        for block in genome:
            if block>=0:
                sequence+=block_library[block+1]
            else:
                sequence+=reverse_complement(block_library[-1*block+1])
    f = open(f"./syntethic/regenerated/{plasmid}_{changes}_{nrOfFragments}.fna","w")
    f.write(sequence)
    f.close()

def write_test(plasmid_1,file,changes,nrOfFragments):
    """Writes the rearranged genome in the "(1 2 ...)" block-string
    format (parsed back by break_distance.readGenomeFromFile()), to
    ./syntethic{file}_{changes}_{nrOfFragments}."""
    sequence = ""

    for genome in plasmid_1:
        sequence+="("
        for block in genome:
            if block>=0:
                sequence+="+"+str(block)+" "
            else:
                sequence+=str(block)+" "
        sequence+=")"

    f = open(f"./syntethic{file}_{changes}_{nrOfFragments}","w")
    f.write(sequence)
    f.close()

def createString(plasmid):
    """Serialises a genome to a plain string, used as a dict/set key in
    generate() to detect when an event produces a genome already seen."""
    pls = ""
    for i in plasmid:
        pls+=str(i)
    return pls 

def generate(plasmid_1,block_library,file,nrOfFragments):
    """Produces k=10 independent rearranged variants of plasmid_1. For each,
    picks a random number of attempted events (0-99) and applies
    REVERSAL/TRANSPOSITION/INSERTION/DELETION events one at a time,
    skipping (not counting) any event whose result exactly reproduces a
    genome already seen in this run — so `changes` undercounts attempted
    events but stays a set of distinct intermediate genomes. The actual
    reported ground-truth distance (changesTrue) is recomputed from
    scratch via TwoBreakDistance between the original and final genome,
    since applied events can partially cancel each other out. Writes one
    (write_test(), write_fasta()) pair per variant."""
    k=10
    events = [Events.REVERSAL,Events.TRANSPOSITION,Events.DELETION,Events.INSERTION]
    curios = []
    for _ in range(k):
        changes = 0
        plasmid = copy.deepcopy(plasmid_1)
        nr_of_changes = random.randrange(100)
        plasmidsCreated = set()
        plasmidsCreated.add(createString(plasmid))
        for j in range(nr_of_changes):
            event = random.choice(events)
            if event == Events.REVERSAL:
                plasmidNew = reversal(plasmid)
            elif event == Events.TRANSPOSITION:
                plasmidNew = transposition(plasmid)
            elif event == Events.INSERTION:
                plasmidNew = insertion(plasmid)
            elif event == Events.DELETION:
                plasmidNew = deletion(plasmid)
            # elif event == Events.TRANSLOCATION:
            #     plasmidNew = translocation(plasmid)
            
            if createString(plasmidNew) not in plasmidsCreated:
                changes+=1
                plasmidsCreated.add(createString(plasmidNew))
                plasmid = plasmidNew

        changesTrue = TwoBreakDistance(plasmid_1,plasmid).dist
        changes = changesTrue
        write_test(plasmid,file,changes,nrOfFragments)
        write_fasta(plasmid,block_library,file,changes,nrOfFragments)

def split(nr,file,nrOfFragments):
    """Splits the single-record FASTA at path "."+file into nrOfFragments
    roughly-equal chunks by line count. Returns [header, block1, block2,
    ...] — this is the block_library that write_fasta() reassembles from,
    and the block IDs 0..nrOfFragments-1 are what run() uses as the
    starting (unrearranged) genome."""
    i = 0
    with open("."+file) as f:
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

def run(file,nrOfFragments):
    """Entry point: splits the plasmid at `file` (path relative to the
    project root, e.g. "/fastas/X.fna") into nrOfFragments blocks and
    generates k=10 synthetic rearranged variants of it via generate()."""
    nr = subprocess.check_output(f"wc -l .{file}",shell=True,text=True)
    nr = nr.strip().split(" ")[0]

    block_library = split(int(nr),file,nrOfFragments)
    plasmid_1 = [[i for i in range(nrOfFragments)]]
    generate(plasmid_1,block_library,file,nrOfFragments)
