import sys
import argparse

def parseArgument():
        # Parse the input
	parser=argparse.ArgumentParser(description=\
                        "Make a script that will convert chromosome names from a bed file to a different naming system for a list of bed files")
	parser.add_argument("--bedFileNameListFileName", required=True,\
                        help='List of bed files with original chromosome names')
	parser.add_argument("--chromNameDictFileNameListFileName", required=True,\
                        help='List of files that map the original chromosome names to the new chromosome names')
	parser.add_argument("--chromNameDictReverse", action="store_true", required=False,\
                        help='Treat the new chomosome names as original chromosome names and the original chromosome names as new chromosome names')
	parser.add_argument("--gzip", action="store_true", required=False,\
                        help='The input and output files are gzipped')
	parser.add_argument("--includeMissingChrom", action="store_true", required=False,\
                        help='Include chromosomes that are not in the dictionary as chrNone')
	parser.add_argument("--outputFileNameSuffix", required=True,\
                        help='Suffix of bed file name where the selected summits will be written, should not start with _, should end with .gz if using gzip')
	parser.add_argument("--numFileNameElementsToRemove", type=int, required=False, default=2,\
			help='Number of parts of the input file name, separted by ., to remove when creating the output file name')
	parser.add_argument("--scriptFileName", required=True,\
                        help='Name of file where the script will be written')
	parser.add_argument("--codePath", required=False, default="/home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode",\
                        help='Path to convertChromNames.py, should not end with /')
	options = parser.parse_args()
	return options

def makeConvertChromNamesScript(options):
	# Make a script that will convert chromosome names from a bed file to a different naming system for a list of bed files
	bedFileNameListFile = open(options.bedFileNameListFileName)
	chromNameDictFileNameListFile = open(options.chromNameDictFileNameListFileName)
	scriptFile = open(options.scriptFileName, 'w+')
	for bedFileNameStr, chromNameDictFileNameStr in zip(bedFileNameListFile, chromNameDictFileNameListFile):
		# Iterate through the bed and corresponding chromosome dictionary files and create a line in the script for each pair
		bedFileName = bedFileNameStr.strip()
		bedFileNameElements = bedFileName.split(".")
		outputFileName =\
			".".join(bedFileNameElements[0:len(bedFileNameElements) - options.numFileNameElementsToRemove]) + "_" + options.outputFileNameSuffix
		scriptFile.write(" ".join(["python", options.codePath + "/convertChromNames.py", "--bedFileName", bedFileName, "--chromNameDictFileName", \
			chromNameDictFileNameStr.strip(), "--outputFileName", outputFileName]))
		if options.chromNameDictReverse:
			# Add the chromNameDictReverse option
			scriptFile.write(" --chromNameDictReverse")
		if options.gzip:
			# Add the gzip option
			scriptFile.write(" --gzip")
		if options.includeMissingChrom:
			# Add the includeMissingChrom option
			scriptFile.write(" --includeMissingChrom")
		scriptFile.write("\n")
	bedFileNameListFile.close()
	chromNameDictFileNameListFile.close()
	scriptFile.close()

if __name__=="__main__":
	options = parseArgument()
	makeConvertChromNamesScript(options)
