import argparse

def parseArgument():
        # Parse the input
        parser = argparse.ArgumentParser(description="Make a script that will run FIMO on a list of fasta files")
        parser.add_argument("--fastaFileNameListFileName", required=False, \
                help="List of fasta files that will be scanned for motifs")
        parser.add_argument("--outputDirSuffix", required=False, default = ".FIMOResults", \
                help="Suffix of directory name where FIMO results will be recorded")
        parser.add_argument("--thresh", required=False, type=float, default=0.001, help="p-Value threshold for FIMO")
        parser.add_argument("--qvThresh", required=False, action="store_true", help="Use a q-Value threshold for FIMO")
	parser.add_argument("--outputToFile", required=False, action="store_true", help="Output the results to a file instead of to a directory")
        parser.add_argument("--scriptFileName", required=True, help="Name of file where script will be written")
	parser.add_argument("--PWMFileName", required=True, help="Name of motif database file")
	
        options = parser.parse_args()
        return options

def makeFIMOListScript(options):
	# Make a script that will run FIMO on a list of fasta files
	fastaFileNameListFile = open(options.fastaFileNameListFileName)
	scriptFile = open(options.scriptFileName, 'w+')
	for line in fastaFileNameListFile:
		# Iterate through the fasta file names and add a line for each to the script
		fastaFileName = line.strip()
		fastaFileNameElements = fastaFileName.split(".")
		fastaFileNamePrefix = ".".join(fastaFileNameElements[0:-1])
		FIMOOutputFileName = fastaFileNamePrefix + options.outputDirSuffix
		scriptFile.write(" ".join(["fimo --max-stored-scores 10000000 --thresh", str(options.thresh)]))
		if options.qvThresh:
			# Add the --qv-thresh option
			scriptFile.write(" --qv-thresh ")
		else:
			# Add the --no-qvalue option
			scriptFile.write(" --no-qvalue ")
		if not options.outputToFile:
			# Add the --o option
			scriptFile.write(" --o " + FIMOOutputFileName + " ")
		else:
			# Add the --text option
			scriptFile.write(" --text ")
		scriptFile.write(" ".join([options.PWMFileName, fastaFileName]))
		if options.outputToFile:
			# Add an output file
			scriptFile.write(" > " + FIMOOutputFileName)
		scriptFile.write("\n")
	fastaFileNameListFile.close()
	scriptFile.close()

if __name__=="__main__":
   options = parseArgument()
   makeFIMOListScript(options)
