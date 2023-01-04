#!/usr/bin/env python3
#foo.py [file with p-vals] [number of folders] [folders with random results] [output file] [max # of trials to consume] [file for extra trials]
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
inFile.readline() 
header = '"Enhancer","Pvalue","Coeff"\n'
new_header = '"Enhancer","Pvalue","Coeff","Exp_Pvalue","Trials","Missing_Trials","Low_Trials"\n'
names = []
pvals = []
parsed_pvals = {}
#corPvals = []
coeffs = []
coeffs_negative = {}
lower_trials = {}
count_trials = {}
for line in inFile:
    tokens = line.strip().replace("\"", "").split(",")
    if len(tokens) >= 2:
        names.append(tokens[0])
        pvals.append(tokens[1])
        parsed_pvals[tokens[0]] = parseNum(tokens[1])
        #corPvals.append(tokens[2])
        coeffs.append(tokens[2])
        coeffs_negative[tokens[0]] = tokens[2][0] == "-"
        if len(tokens) >= 7:
            lower_trials[tokens[0]] = int(tokens[6])
            count_trials[tokens[0]] = int(tokens[4])
        else:
            lower_trials[tokens[0]] = 1
            count_trials[tokens[0]] = 1
            
inFile.close()

#lower_trials = {}
#count_trials = {}
#for name in names:
#    lower_trials[name] = 1
#    count_trials[name] = 1

nDir = int(sys.argv[2])

max_trials = int(sys.argv[4+nDir])
extraFile = open(sys.argv[5+nDir], "w")
extraFile.write(header)
 
for i in range(nDir):   
    inDir = sys.argv[3+i]
    print(inDir)
    for f in os.listdir(inDir):
        inFile = open(inDir + "/" + f, "r") 
        header = inFile.readline().strip().replace("\"", "").split(",")
        try:
            column = header.index("Pvalue")
            column2 = header.index("Coeff")
        except:
            print(f, file=sys.stderr)
            continue
        for line in inFile:
            tokens = line.strip().replace("\"", "").split(",")
            name = tokens[0]
            if name not in coeffs_negative:  #Skip for case where perms contain unfiltered peaks
                continue
            if count_trials[name] == max_trials:
                extraFile.write(line)
                continue
            if (tokens[column2] != "0" and (tokens[column2][0] == "-") != coeffs_negative[name]) or (tokens[column2] == "0" and parsed_pvals[name][0] != 0):
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
    missing_trials = 10 ** max(3, ceil(log10(trials))) - trials
    low_trials = lower_trials[name]
    outFile.write("," + str(lower_trials[name] / count_trials[name]) + "," + str(trials) + "," + str(missing_trials) + "," + str(low_trials) + "\n")
outFile.close()
extraFile.close()
