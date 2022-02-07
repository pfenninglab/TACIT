#!/usr/bin/env python3

#python3 p_dist_.py gene-peak-file gene-p-val-file outFile


import sys
from bisect import insort
from scipy.stats import ranksums




peaksFile = open(sys.argv[1], "r")
fgGenes = set()
for line in peaksFile:
    tokens = line.strip().split()
    peak = tokens[3]
    gene = tokens[7]
    fgGenes.add(gene)
peaksFile.close()

fgPvals = []
bgPvals = []
pValFile = open(sys.argv[2], "r")
pValFile.readline()
for line in pValFile:
    tokens = line.strip().split("\t")
    gene = tokens[0].replace("\"", "")
    if tokens[3] == "NA":
        continue
    p = float(tokens[3])
    if gene in fgGenes:
        insort(fgPvals, p)
    else:
        insort(bgPvals, p)
pValFile.close()

print(len(fgPvals), len(bgPvals))

print(ranksums(fgPvals, bgPvals, "less"))

outFile = open(sys.argv[3], "w")
outFile.write("pval,fg\n")
for pval in fgPvals:
    outFile.write(str(pval) + ",T\n")
for pval in bgPvals:
    outFile.write(str(pval) + ",F\n")
outFile.close()

