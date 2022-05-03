import sys
import argparse
import gzip

def parseArgument():
	# Parse the input
	parser=argparse.ArgumentParser(description=\
			"Convert chromosome names from a bed file to a different naming system")
	parser.add_argument("--bedFileName", required=True,\
			help='Bed file with original chromosome names')
	parser.add_argument("--chromNameDictFileName", required=True,\
                        help='File that maps the original chromosome names to the new chromosome names')
	parser.add_argument("--chromNameDictReverse", action="store_true", \
			required=False,\
                        help='Treat the new chomosome names as original chromosome names and the original \
				chromosome names as new chromosome names')
	parser.add_argument("--removeSuffixes", action="store_true", \
                        required=False,\
                        help='Remove parts after . in bed file chromosome names')
	parser.add_argument("--gzip", action="store_true", required=False,\
                        help='The input and output files are gzipped')
	parser.add_argument("--includeMissingChrom", action="store_true", \
			required=False,\
                        help='Include chromosomes that are not in the dictionary as NA')
	parser.add_argument("--outputFileName", required=True,\
			help='bed file where the regions with converted chromosome names will be written')
	options = parser.parse_args()
	return options

def convertChromNames(options):
	# Convert chromosome names from a bed file to a different naming system
	chromNameDict = {}
	chromNameDictFile = open(options.chromNameDictFileName)
	for line in chromNameDictFile:
		# Iterate through the chromosomes and create an entry in the dictionary for each
		lineElements = line.strip().split("\t")
		if lineElements[1] == "-":
			# There is no chromosome in the conversion
			continue
		if options.chromNameDictReverse:
			# Make the second column the keys and the first column the values
			chromNameDict[lineElements[1]] = lineElements[0]
		else:
			# Make the first column the keys and the second column the values
			chromNameDict[lineElements[0]] = lineElements[1]
	chromNameDictFile.close()
	bedFile =\
		gzip.open(options.bedFileName) if options.gzip else \
			open(options.bedFileName) # Use gzip to open the bed file if the bed file is gzipped
	outputFile =\
		gzip.open(options.outputFileName, 'wb') if options.gzip else \
			open(options.outputFileName, 'w+') # Use gzip to open the output file if the output file should be gzipped
	for line in bedFile:
		# Iterate through the lines of the peak motif hit file and select peak summits with a near-by motif hit
		lineElements = line.strip().split("\t")
		chrom = lineElements[0]
		if options.removeSuffixes:
			# Remove the part of the chromosome name after the .
			chromElements = chrom.split(".")
			chrom = chromElements[0]
		if chrom not in chromNameDict:
			# Skip this region or make it chrNone because the chromosome is not in the list of chromosomes
			if options.includeMissingChrom:
				# Include the region as NA
				outputFile.write("NA" + "\t" + \
					"\t".join(lineElements[1:len(lineElements)]) + "\n")
				continue
			print("Problem: Chromosome " + chrom + \
				" not in list of chromosomes")
			continue
		outputFile.write(chromNameDict[chrom] + "\t" + \
			"\t".join(lineElements[1:len(lineElements)]) + \
			"\n")
	bedFile.close()
	outputFile.close()

if __name__=="__main__":
	options = parseArgument()
	convertChromNames(options)
