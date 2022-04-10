#Collects specific p-values for OCR-gene pairs. Takes four arguments:
#1) File with OCR coordinates chr start end OCR
#2) File with OCR, Gene1, Gene2, Gene3, ...
#3) File with TSSs, chr \t TSS \t __ \ Gene1
#4) HiC p-values file (gzip, p-value in 6th column)
#Output: tab-seperated OCR Gene Pvalue #Contacts Note

import sys, gzip

#all OCRs -> [chr, start, end]
ocrDict = {}
coordFile = open(sys.argv[1], "r")
for line in coordFile:
    tokens = line.split()
    ocrDict[tokens[3]] = tokens[:3]

geneFile = open(sys.argv[2], "r")
#specific OCRs -> [genes]
geneDict = {}
ocrList = [] #keep order
for line in geneFile:
    tokens = [x.strip() for x in line.split(",")]
    geneDict[tokens[0]] = tokens[1:]
    ocrList.append(tokens[0]) 
geneFile.close()


#chr_bin -> OCR
binDict = {}
for ocr in geneDict:
    if ocr in ocrDict:
        chr, start, end = ocrDict[ocr]
        middle = (int(start) + int(end)) // 2
        bin = str(chr) + "_" + str(middle)[:-4] + "5000" #Hardcoded
        if bin in binDict:
            print("Two OCRs in bin:", binDict[bin], ocr)
        else:
            binDict[bin] = ocr

#gene -> TSS bins (chr_bin)
tssDict = {}
tssFile = open(sys.argv[3], "r")
for line in tssFile:
    tokens = line.split()
    gene = tokens[3]
    bin = tokens[0] + "_" + tokens[1][:-4] + "5000"
    if gene not in tssDict:
        tssDict[gene] = set()
    tssDict[gene].add(bin)


#ocr -> {gene -> pval}
pvalDict = {}
with gzip.open(sys.argv[4], 'rb') as f:
    file_content = f.read()
hicLines = file_content.decode().split("\n")[1:]
for line in hicLines:
    tokens = line.strip().split()
    if len(tokens) == 0:
        continue
    bins = [tokens[0] + "_" + tokens[1], tokens[2] + "_" + tokens[3]]
    p = float(tokens[5])
    cont = int(tokens[4])
    for i in range(2):
        bin = bins[i]
        otherBin = bins[(i+1) %2]
        if bin in binDict:
            ocr = binDict[bin]
            if ocr not in pvalDict:
                pvalDict[ocr] = {}
            genes = geneDict[ocr]
            for gene in genes:
                if gene in tssDict:
                    tss = tssDict[gene]
                    if otherBin in tss:
                        if gene not in pvalDict[ocr] or p < pvalDict[ocr][gene][0]:
                            pvalDict[ocr][gene] = (p, cont)

print("OCR\tGene\tPvalue\tContacts\tNote")
for ocr in ocrList:
    for gene in geneDict[ocr]: 
        note = ""
        cont = 0
        if ocr in pvalDict:
            if gene in pvalDict[ocr]:
                p, cont = pvalDict[ocr][gene]
                
            else:
                p = 1
                note = "Gene Not Found"
        else:
            p = -1
            note = "No OCR Ortholog"
        print(ocr, gene, p, cont, note, sep="\t")
      


    