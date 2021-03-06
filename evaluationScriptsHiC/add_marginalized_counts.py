#Adds marginalized counts to a list of fragments (bins):
#1) fragments file generated by createFitHiCFragments-fixedsize.py (gzipped)
#2) contacts file generated by validPairs2FitHiC-fixedSize.sh (gzipped)
# with default values (1) in 4th and fifth column
#Output: fragments file, with fourth column marginalized counts and
#Fifth column 1 if marginalized counts >= 0 and 0 otherwise

import sys, gzip

fragment_list = []
count_dict = {}

with gzip.open(sys.argv[1], 'rb') as f:
    file_content = f.read()
lines = file_content.decode().split("\n")
for line in lines:
    tokens = line.strip().split()
    if len(tokens) != 5:
        continue
    fragment_list.append(tokens[:3])
    count_dict[tokens[0] + "_" + tokens[2]] = 0

with gzip.open(sys.argv[2], 'rb') as f:
    file_content = f.read()
lines = file_content.decode().split("\n")
for line in lines:
    tokens = line.strip().split()
    if len(tokens) != 5:
        continue
    fragment1 = tokens[0] + "_" + tokens[1]
    fragment2 = tokens[2] + "_" + tokens[3]
    pairCount = int(tokens[4])
    count_dict[fragment1] += pairCount
    count_dict[fragment2] += pairCount

for f in fragment_list:
    count = count_dict[f[0] + "_" + f[2]]
    f.append(str(count))
    f.append("1" if count > 0 else "0")
    print("\t".join(f))

    




 