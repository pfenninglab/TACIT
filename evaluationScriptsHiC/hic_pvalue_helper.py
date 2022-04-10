#Assists in manual collection of HiC p-values from Hugin. Takes four-five arguments:
#1) File with OCR coordinates (chr start end OCR)
#2) File with (OCR, Gene1, Gene2, Gene3, ...)
#3) File with TSSs (chr \t TSS \t __ \ Gene1)
#4) File to write collected data:
#   (OCR chr start end TSS start end count exp, -log10(p)
#5) Restart file - file containing some previously-collected data,
#   i.e. argument 4 of a previous run
#Output written to STDOUT

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


#gene -> TSS bins (chr_bin)
tssDict = {}
tssFile = open(sys.argv[3], "r")
for line in tssFile:
    tokens = line.split()
    gene = tokens[3]
    bin = tokens[1][:-4] + "0000"
    if gene not in tssDict:
        tssDict[gene] = set()
    tssDict[gene].add(bin)
tssFile.close()
    
partialDict = {}
if len(sys.argv) > 5:
    partialFile = open(sys.argv[5])
    partialFile.readline()
    for line in partialFile:
        tokens = line.strip().split()
        if len(tokens) == 8:
            key = tokens[0] + "_" + tokens[1] + "_" + tokens[3]
            partialDict[key] = tokens[5:8]
    partialFile.close()
    
#ocr -> {gene -> pval}
pvalDict = {}

outFile = open(sys.argv[4], "w")
outFile.write("\t".join(("Chr", "Start", "End", "Start", "End", "Contacts", "Expected", "-log10(p)")) + "\n")

for ocr in ocrList:
    if ocr not in ocrDict:
        continue
    print("\n\nNext OCR:", ocr)
    if ocr in []: #List of OCRs known to have no contacts
        pvalDict[ocr] = "EMPTY"
        continue
    if ocr in []: #List of OCRs known to have too few contacts
        pvalDict[ocr] = "FEW"
        continue
    chr, start, end = ocrDict[ocr]
    if chr == "chrX":
        pvalDict[ocr] = "SEX"
        continue
    print("Anchor:")
    middle = (int(start) + int(end)) // 2
    print(chr + ":" + str(middle))
    print("Window:")
    print(chr+":"+str(middle-2000000)+"-"+str(middle+2000000))
    ocr_start = str(middle)[:-4] + "0000"
    ocr_end = str(middle+10000)[:-4] + "0000"
    pvalDict[ocr] = {}
    genes = geneDict[ocr]
    for gene in genes:
        if gene in tssDict:
            tsss = tssDict[gene]
            for tss in tsss:
                partialKey = chr + "_" + ocr_start + "_" + tss
                if partialKey in partialDict:
                    cont, exp, p = partialDict[partialKey]
                else:
                    print("\nOther bin:", tss)
                    cont = input("Contacts: ")
                    exp = input("Expected: ")
                    p = input("-log10 p: ")
                outFile.write("\t".join((chr, ocr_start, ocr_end, tss, str(int(tss)+10000), cont, exp, p)) + "\n")
                cont = int(cont)
                p = 10 ** (-float(p))
                if gene not in pvalDict[ocr] or p < pvalDict[ocr][gene][0]:
                    pvalDict[ocr][gene] = (p, cont)
    

outFile.close()

print("OCR\tGene\tPvalue\tContacts\tNote")
for ocr in ocrList:
    for gene in geneDict[ocr]: 
        note = ""
        cont = 0
        if ocr in pvalDict:
            if pvalDict[ocr] == "SEX":
                p = -1
                note = "OCR Ortholog on sex chromosome"            
            elif pvalDict[ocr] == "FEW":
                p = -1
                note = "OCR excluded bc few nearby contacts"
            elif pvalDict[ocr] == "EMPTY":
                p = -1
                note = "OCR excluded bc no nearby contacts"
            elif gene in pvalDict[ocr]:
                p, cont = pvalDict[ocr][gene]
                
            else:
                p = 1
                note = "Gene Not Found"
        else:
            p = -1
            note = "No OCR Ortholog"
        print(ocr, gene, p, cont, note, sep="\t")
outFile.close()
      


    