#!/usr/bin/env python3

#python3 phylop_p_correlation.py [OCR BED file], [p-value file], [human | mouse]
# Human & mouse PhyloP paths are hardcoded. Assumes per-chromosome file for human
# and a single file for mouse. OCRs in BED file must be a subset of those in p-value file.

import pyBigWig, sys
from scipy.stats import pearsonr

HUMPATH = "/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/phyloP/human-centered-200m-Feb2021/"
MOUSEPATH = "/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseGenome/"

#Parse BED file
coordDict = {}
bedFile = open(sys.argv[1])
for line in bedFile:
    tokens = line.split()
    coordDict[tokens[3]] = [tokens[0], int(tokens[1]), int(tokens[2])]
bedFile.close()

#Parse OCR-phenotype p-values or similar
pDict = {}

resFile = open(sys.argv[2])
resFile.readline()
for line in resFile:
    tokens = line.split(",")
    pDict[tokens[0]] = float(tokens[1])
resFile.close()

#Obtain average PhyloPs
pvals = []
phylops = []
if sys.argv[3] == "human":
    chrDict = {}
    chrs = [str(x) for x in range(1,23)] + ["X", "Y"]
    for chr in chrs:
        chrDict["chr" + chr] = []
    for ocr in pDict:
        if ocr in coordDict:
            chr = coordDict[ocr][0]
            if chr in chrDict:
                chrDict[chr].append(ocr)
    for chr in chrs:
        print(chr)
        bwFile = pyBigWig.open(HUMPATH + f"200m_scoresPhyloP_20210214.{chr}.bigWig")
        chrName = list(bwFile.chroms())[0]
        for ocr in chrDict["chr" + chr]:
            start, end = coordDict[ocr][1:]
            pvals.append(pDict[ocr])
            phylops.append(bwFile.stats(chrName, start, end, exact=True)[0])
            
elif sys.argv[3] == "mouse":
    i = 0
    bwFile = pyBigWig.open(MOUSEPATH + "phylop.bw")
    chrNames = bwFile.chroms().keys()
    for ocr in pDict:
        if ocr in coordDict:
            chrName, start, end = coordDict[ocr]
            if chrName in chrNames:
                pvals.append(pDict[ocr])
                phylops.append(bwFile.stats(chrName, start, end, exact=True)[0])
            else:
                print(chrName)
        i += 1
        if (i % 1000) == 0:
            print(i)

print("\n\nCorrelation:", pearsonr(pvals, phylops))


    


    

