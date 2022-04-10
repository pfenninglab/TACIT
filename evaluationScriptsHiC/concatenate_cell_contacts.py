#Concatenates gzipped contacts from single cells
#Argument: directory containing all and only pairs.txt.gz files

import os, sys, gzip

files = os.listdir(sys.argv[1])
first = True
for name in files:
    with gzip.open(sys.argv[1] + "/" + name, 'rb') as f:
        file_content = f.read()
    lines = file_content.decode().split("\n")
    for line in lines:
        if len(line) == 0:
            print()
            continue
        if line[0] == "#":
            if first:
                print(line)
        else:
            print(line)   
    first = False
    