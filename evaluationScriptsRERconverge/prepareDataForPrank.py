import sys
import argparse
import os
import gzip
import pyfastx

def parseArgument():
        # Parse the input
        parser=argparse.ArgumentParser(description=\
                        "Prepare data for running PRANK")
        parser.add_argument("--bedFileNameListFileName", required=True,\
                        help='List of bed files from different species, should have peak names and, if extendedSummits option is being used, be in narrowPeak format')
        parser.add_argument("--baseSpeciesLine", type=int, required=True,\
                        help='Line in the species file with the base species, should be 0-indexed')
        parser.add_argument("--peaksToIncludeFileName", required=False,\
                        help='Name of file with list of names of peaks to include, will include all peaks if this is not provided')
        parser.add_argument("--fastaFileNameListFileName", required=False,\
                        help='Name of file with list of fasta files from each species, should be in the same order as bedFileNameListFileName')
        parser.add_argument("--chromSizesFileNameListFileName", required=False,\
                        help='Name of file with list of chromosome sizes each species, should be in the same order as bedFileNameListFileName')
        parser.add_argument("--extendedSummits", action="store_true", required=False,\
                        help='Fasta sequences are from extended summits instead of peaks, assumes that fastaFromPeaks option is being used')
        parser.add_argument("--sequenceLength", type=int, required=False, default=500,\
                        help='Total sequence length sequences are summit extensions, summits will be extended sequenceLength/2, assumes extendedSummits is being used')
        parser.add_argument("--speciesListFileName", required=True,\
                        help='Name of file with species list')
        parser.add_argument("--speciesTreeFileName", required=True,\
                        help='Name of file with species tree')
        parser.add_argument("--prankDataDir", required=True,\
                        help='Directory where inidividual fasta files for each peak and prank outputs will be stored')
        parser.add_argument("--individualPeakFileNameSuffix", required=False, default="orthologs.fa", \
                        help='Suffix that will be added to fasta ortholog files in addition to base species and its peak name and peak region')
        parser.add_argument("--outputFileNameSuffix", required=False, default="prankOutputs.fas", \
                        help='Suffix that will be added to prank outputs in addition to base species and its peak name and peak region')
        parser.add_argument("--scriptFileName", required=True, \
                        help='File where all the prank commands will be recorded')
        options = parser.parse_args()
        return options

def prepareDataForPrank(options):
	bedFileNameListFile = open(options.bedFileNameListFileName)
	chromSizesFileNameListFile = open(options.chromSizesFileNameListFileName)
	bedDictList = []
	extensionSize = int(options.sequenceLength/2)
	peaksToInclude = []
	if options.peaksToIncludeFileName != None:
                # Only a strict subset of peaks will be included
                peaksToIncludeFile = open(options.peaksToIncludeFileName)
                for line in peaksToIncludeFile:
                        # Iterate through the list of peaks to include and add each to the list of peaks to include
                        peaksToInclude.append(line.strip())
                peaksToIncludeFile.close()
                print ("Inclusion peak list made!")
	for bedFileNameStr, chromSizesFileName in zip(bedFileNameListFile, chromSizesFileNameListFile):
		# Iterate through the bed files and make a dictionary mapping peaks to coordinates for each
		bedFileName = bedFileNameStr.strip()
		bedFile = None
		if bedFileName.endswith("gz"):
			# The bed file is gzipped, so use gzip to open it
			bedFile = gzip.open(bedFileName, 'rt')
		else:
			# Open the bed file without gzip
			bedFile = open(bedFileName, 'rt')
		chromSizesFile = open(chromSizesFileName.strip())
		chromSizesDict = {}
		for line in chromSizesFile:
			# Iterate through the chromosome sizes and make a dictionary entry for each chromosome
			lineElements = line.strip().split("\t")
			chromSizesDict[lineElements[0]] = int(lineElements[1])
		chromSizesFile.close()
		bedDict = {}
		for line in bedFile:
			# Iterate through the lines of the bed file and make an entry in the dictionary for each
			lineElements = line.strip().split("\t")
			chrom = lineElements[0]
			if chrom not in chromSizesDict:
				# Do not use the current chromosome because its size is unknown
				print("Chromosome " + chrom + " is not in " + chromSizesFileName.strip())
				continue
			peakName = lineElements[3]
			if (len(peaksToInclude) > 0) and (peakName not in peaksToInclude):
				# The current peak should not be included
				continue
			start = None
			end = None
			if options.extendedSummits:
				# Use region surrounding summits instead of given region
				start = int(lineElements[1]) + int(lineElements[9]) - extensionSize + 1
				end = int(lineElements[1]) + int(lineElements[9]) + extensionSize
			else:
				# Use the peak
				start = int(lineElements[1])
				end = int(lineElements[2])
			if (start < 0) or (end > chromSizesDict[chrom]):
				# The current region is outside of the chromosome's bound
				continue
			bedDict[peakName] = (chrom, start, end)
		bedFile.close()
		bedDictList.append(bedDict)
	print("Number of bed files is " + str(len(bedDictList)))
	print ("Peak name to coordinate dictionaries made!")
	bedFileNameListFile.close()
	chromSizesFileNameListFile.close()
	fastaFileNameListFile = open(options.fastaFileNameListFileName, 'rt')
	fastaList = []
	fastaChromList = []
	for line in fastaFileNameListFile:
		# Iterate through the fasta files and add each to a list
		fastaFileName = line.strip()
		assert(os.path.exists(fastaFileName))
		fa = pyfastx.Fasta(fastaFileName)
		fastaChroms = set(fa.keys())
		fastaList.append(fa)
		fastaChromList.append(fastaChroms)
	fastaFileNameListFile.close()
	speciesListFile = open(options.speciesListFileName, 'rt')
	speciesList = [line.strip() for line in speciesListFile]
	print("Number of species is " + str(len(speciesList)))
	numSpecies = len(speciesList)
	speciesListFile.close()
	scriptFile = open(options.scriptFileName, 'w+')
	print ("Species name and assembly lists made!")
	for peak in bedDictList[options.baseSpeciesLine].keys():
                # Iterate through the peaks and prepare each one for PRANK
                peakInfoStr =\
			speciesList[options.baseSpeciesLine] + "_" + \
				"_".join([bedDictList[options.baseSpeciesLine][peak][0], str(bedDictList[options.baseSpeciesLine][peak][1]), \
					str(bedDictList[options.baseSpeciesLine][peak][2]), peak])
                individualPeakFastaFileName = options.prankDataDir + "/" + peakInfoStr + "_" + options.individualPeakFileNameSuffix
                individualPeakFastaFile = open(individualPeakFastaFileName, 'w+')
                for i in range(numSpecies):
			# Iterate through the species and get the sequences of the current individual peak in each
                        bedDict = bedDictList[i]
                        species = speciesList[i]
                        if peak not in bedDict.keys():
				# The current peak does not have a usuable ortholog in the current species
                                if i == numSpecies - 1:
					# At the end of the species list, so stop
                                        break
                                continue
                        fa = fastaList[i]
                        if bedDict[peak][0] not in fastaChromList[i]:
				# The current chromosome is not in the fasta file
                                if i == numSpecies - 1:
					# At the end of the species list, so stop
                                        break
                                continue
                        currentSequence = fa.fetch(bedDict[peak][0], (bedDict[peak][1], bedDict[peak][2]))
                        individualPeakFastaFile.write(">" + species + "\n")
                        individualPeakFastaFile.write(currentSequence + "\n")
                outputFileName = options.prankDataDir + "/" + peakInfoStr + "_" + options.outputFileNameSuffix
                scriptFile.write(" ".join(["prank -d=" + individualPeakFastaFileName, "-t=" + options.speciesTreeFileName, "-o=" + outputFileName, \
			"-F -showxml -showanc"]) + "\n")
                individualPeakFastaFile.close()
	scriptFile.close()

if __name__=="__main__":
        options = parseArgument()
        prepareDataForPrank(options)
