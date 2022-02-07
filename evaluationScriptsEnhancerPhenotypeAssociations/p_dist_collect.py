#!/usr/bin/env python3

#python3 p_dist_.py gene-list gene-peak-file p-val-file outFile


import sys
from bisect import insort
from scipy.stats import ranksums


geneFile = open(sys.argv[1], "r")
geneSet = set([x.strip() for x in geneFile])
geneFile.close()

peaksFile = open(sys.argv[2], "r")
fgPeaks = set()
allPeaks = set()
for line in peaksFile:
    tokens = line.strip().split()
    peak = tokens[3]
    gene = tokens[7]
    allPeaks.add(peak)
    if gene in geneSet:
        fgPeaks.add(peak)
peaksFile.close()

bgPeaks = allPeaks - fgPeaks

fgPvals = []
bgPvals = []
pValFile = open(sys.argv[3], "r")
pValFile.readline()
for line in pValFile:
    tokens = line.strip().split(",")
    p = float(tokens[3])
    if tokens[0] in fgPeaks:
        insort(fgPvals, p)
    elif tokens[0] in bgPeaks:
        insort(bgPvals, p)
pValFile.close()

print(len(fgPvals), len(bgPvals))

print(ranksums(fgPvals, bgPvals, "less"))

outFile = open(sys.argv[4], "w")
outFile.write("pval,fg\n")
for pval in fgPvals:
    outFile.write(str(pval) + ",T\n")
for pval in bgPvals:
    outFile.write(str(pval) + ",F\n")
outFile.close()

