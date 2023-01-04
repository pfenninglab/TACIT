#!/usr/bin/env python3

#foo.py [file with p-vals] [number of folders] [folders with random results] [output file]

#Script to compute permulation p-values from ocr_phylo(g)lm.r results,
#using the method of Kulinskaya, 2008, in which permulations in the other
#direction than the initial trial are rejected.
#Takes at least four required positional arguments in order:
#   The result CSV file of ocr_phylo(g)lm.r run once without permulations for every OCR,
#   The number of directories containing results files
#   Directories containing results file(s) of ocr_phylo(g)lm run with permulations on
#       a subset of the OCRs in the first file. The number of directories provided should
#       equal the previous argument.
#   An output CSV file, which will be a copy of the input file with additional columns giving
#       the permulations p-value and number of trials performed for each OCR


import os, sys
from math import log10, ceil

def parseNum(s):
    te = s.split("e")
    if len(te) == 2:
        num = float(te[0])
        base = int(te[1])
    else:
        num = float(s)
        base = 0
    return [num, base]

def leq(a, b):
    randNum, randBase = a
    realNum, realBase = b
    return randBase < realBase or (realBase == randBase and randNum <= realNum)


inFile = open(sys.argv[1], "r")

new_header = inFile.readline().strip() + ",\"Exp_Pvalue\",\"Trials\",\"Missing_Trials\"\n" #Name, P-val, Coeff, Exp P-val, Trials
names = []
pvals = []
parsed_pvals = {}
#corPvals = []
coeffs = []
coeffs_negative = {}
for line in inFile:
    tokens = line.strip().replace("\"", "").split(",")
    if len(tokens) == 3:
        names.append(tokens[0])
        pvals.append(tokens[1])
        parsed_pvals[tokens[0]] = parseNum(tokens[1])
        #corPvals.append(tokens[2])
        coeffs.append(tokens[2])
        coeffs_negative[tokens[0]] = tokens[2][0] == "-"
inFile.close()

lower_trials = {}
count_trials = {}
for name in names:
    lower_trials[name] = 1
    count_trials[name] = 1

nDir = int(sys.argv[2])   
for i in range(nDir):   
    inDir = sys.argv[3+i]
    for f in os.listdir(inDir):
        inFile = open(inDir + "/" + f, "r") 
        header = inFile.readline().strip().replace("\"", "").split(",")
        column = header.index("Pvalue")
        column2 = header.index("Coeff")
        #print(f)
        for line in inFile:
            tokens = line.strip().replace("\"", "").split(",")
            name = tokens[0]
            if ((tokens[column2][0] == "-") != coeffs_negative[name]) or tokens[column2] == "0":
                continue #Reject
            pval = tokens[column]
            count_trials[name] = 1 + count_trials[name]
            if leq(parseNum(pval), parsed_pvals[name]):
                lower_trials[name] = lower_trials[name] + 1
        

outFile = open(sys.argv[3+nDir], "w")
outFile.write(new_header)
for i in range(len(names)):
    name = names[i]
    outFile.write(",".join([name, pvals[i], coeffs[i]]))
    trials = count_trials[name]
    missing_trials = 10 ** max(4, ceil(log10(trials))) - trials
    outFile.write("," + str(lower_trials[name] / count_trials[name]) + "," + str(trials) + "," + str(missing_trials) + "\n")
outFile.close()
