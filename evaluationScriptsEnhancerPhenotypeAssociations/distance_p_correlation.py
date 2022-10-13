#!/usr/bin/env python3

#python3 distance_p_correlation.py [OCR BED file], [p-value file]
# OCR BED file should contain distances to the nearest TSS as the last column.
# Test will be performed on OCRs in BED file, which must be a subset of those
# in the p-value file. 

import sys
from math import log10
from scipy.stats import pearsonr


#Parse distance file
distDict = {}
bedFile = open(sys.argv[1])
for line in bedFile:
    tokens = line.split()
    peak = tokens[3]
    dist = int(tokens[-1])
    distDict[peak] = dist
bedFile.close()

#Parse OCR-phenotype p-values or similar
pDict = {}

resFile = open(sys.argv[2])
resFile.readline()
for line in resFile:
    tokens = line.split(",")
    pDict[tokens[0]] = float(tokens[1])
resFile.close()

peaks = list(distDict.keys())
dists = [distDict[peak] for peak in peaks]
#logDists = [log10(x) for x in dists]
pvals = [pDict[peak] for peak in peaks]

print("\n\nCorrelation:", pearsonr(pvals, dists))


    


    

