

########################################################################################
##  This is the script for running an ordinary EM-like algorithm on bipartite graphs  ##
########################################################################################



import operator
from random import randint
import networkx as nx
from math import sqrt, log
from copy import copy
import sys




filename = sys.argv[1]
weights_file = sys.argv[2]
outputfile = sys.argv[3]





reads = {} # reads = {read1: [protein1, protein2, ...]} hashmap storing the list of proteins to which a read maps
proteins = {} # proteins = {protein1: [read1, read2, ...]} hashmap storing the list of reads that mapped to it
prexp = {} # prexp = {protein: exp} "expression" of proteins
new_prexp = {} # this will store the same thing as prexp, but updated after each iteration
names_map = {}
weights_map = {}
#prot_len = {}

INF = 10000000



with open(weights_file) as f:
    weights = f.readlines()
    weights = map(lambda x: x.strip().split(), weights)
    for x, y in weights:
        weights_map[x] = float(y)



def distance(new_prexp, prexp):
    dist = 0
    for key in prexp:
        dist += (prexp[key] - new_prexp[key]) ** 2
    return sqrt(dist)
    





with open(filename) as f:
    for line in f.readlines():
        line = line.split()
        read, protein = line[0][:17], line[2]
        e_value = float(line[9])
        read_where_maps = reads.get(read, set())
        read_where_maps.add((protein, e_value))
        reads[read] = read_where_maps

        protein_which_reads = proteins.get(protein, set())
        protein_which_reads.add(read)
        proteins[protein] = protein_which_reads




for protein in proteins:
    # initialize protein "expressions" by random numbers
    new_prexp[protein] = randint(1, 1000)
    prexp[protein] = INF


while abs(distance(new_prexp, prexp)) > 0.0001:
    prexp = copy(new_prexp)
    # E-step:
    # Compute the expected number n(j) of reads that come from protein j
    nj = {}
    for protein in proteins:
        nj[protein] = 0
    for read in reads:
        s = 0
        for protein, e_value in reads[read]:
            freq_prot = weights_map[read]
            weight = float(e_value) #/ float(e_value) #prot_len[protein] / float(e_value)
            s += freq_prot * weight * prexp[protein]
        for protein, e_value in reads[read]:
            freq_prot = weights_map[read]
            weight = float(e_value)
            if s > 0:
                nj[protein] += freq_prot * weight * prexp[protein] * 1.0 / s
            else:
                nj[protein] = 0
    # M-step:
    suma = 0
    for protein in nj:
        suma += nj[protein]# * 1.0 / prot_len[protein]
    for protein in nj:
        new_prexp[protein] = nj[protein] * 1.0 / suma #(nj[protein] * 1.0 / prot_len[protein]) / suma

    print distance(new_prexp, prexp), "this is the distance"



with open(outputfile, "w") as q:
    filtered = {}
    for x, y in new_prexp.iteritems():
        #if float(y) > 0.000001:
        #    filtered[x] = y
        filtered[x] = y
    filtered = sorted(filtered.items(), key=operator.itemgetter(1))
    for x, y in filtered:
        q.write("%s    %s\n" % (x, y * 10000))



