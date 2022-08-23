#python3 printCommandsFor_getfasta.py <fastaKeyFile.txt> <listOfBedFiles.txt> [numCharToRemoveFront] [numCharToRemoveBack] <OUT_scriptFileName.sh>

import csv
import sys


fastaKeyFile = sys.argv[1]
listFile = sys.argv[2]
r1 = int(sys.argv[3])
r2 = int(sys.argv[4])
scriptFileName = sys.argv[5]

fastaKey = {}
with open(fastaKeyFile, 'r') as fk:
    reader = csv.reader(fk, dialect='excel', delimiter='\t')
    for line in reader:
        fastaKey[line[0]] = line[1]

scriptFile = open(scriptFileName, 'w+')
with open(listFile, 'r') as f:
    reader = csv.reader(f, dialect='excel', delimiter=' ')
    for row in reader:
        speciesName=row[0][r1:-r2]
        if speciesName in fastaKey:
            scriptFile.write("bedtools getfasta -name -fi "+fastaKey[speciesName]+" -bed "+row[0]+" -fo fasta/"+row[0][r1:-4]+".fa"+"\n")
        else:
            print("ERROR: no fasta file for "+speciesName)
