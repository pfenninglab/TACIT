#python3 printCommandsFor_HALPER2bed.py <listOfHALPERoutFiles.txt> [numCharToRemoveFront] [numCharToRemoveBack] [fileSuffix] <OUT_scriptFileName.sh>

import csv
import sys


listFile = sys.argv[1]
r1 = int(sys.argv[2])
r2 = int(sys.argv[3])
suffix = sys.argv[4]
scriptFileName = sys.argv[5]

scriptFile = open(scriptFileName, 'w+')
with open(listFile, 'r') as f:
    reader = csv.reader(f, dialect='excel', delimiter=' ')
    for row in reader:
        speciesName=row[0][r1:-r2]
        scriptFile.write("awk 'OFS=\"\\t\" {print $1,$4-250,$4+250,$5\"-"+speciesName+"\"}' "+row[0]+" > tmp.bed"+"\n")
        scriptFile.write("awk '($2 > 0)' tmp.bed > bedOut/"+speciesName+suffix+".bed"+"\n")
