import sys
import argparse

"""
This script makes a script that will make summit-centered fasta files for list of narrowPeak files and corresponding fastas and chromosome sizes files.
To run: python makeSummitCenteredFastaFromNarrowPeakScript.py [options]
"""

def parseArgument():
	# Parse the input
	parser = argparse.ArgumentParser(description =\
		"Make a script that will make summit-centered fasta files for \
			a list of narrowPeak files")
	parser.add_argument("--narrowPeakFileNameListFileName", \
		required=True, \
		help='File with list of narrowPeak file names')
	parser.add_argument("--maxPeakLength", type=int, required=False, \
		default=None, help='Maximum number of bases per peak')
	parser.add_argument("--sequenceLength", type=int, required=False, \
		default=500, help='Number of bases in each sequence')
	parser.add_argument("--genomeFileNameListFileName", required=True, \
		help='File with genome sequence file list, \
			where the genome in each line corresponds to \
			the narrowPeak file in the same line')
	parser.add_argument("--chromSizesFileNameListFileName", required=True, \
                help='File with chromosome sizes file list, \
			where the chromosome sizes file \
			in each line corresponds to \
			the narrowPeak file in the same line')
	parser.add_argument("--createOptimalBed", action='store_true', \
		required=False, \
                help='Remove peaks on unknown, alternative, \
			and random chromosomes before making fasta files')
	parser.add_argument("--nameOnly", action='store_true', required=False, \
		help = 'Use the nameOnly option when making the fasta file')
	parser.add_argument("--path", required=True, \
		help='Path to makeSummitCenteredFastaFromNarrowPeak.py')
	parser.add_argument("--scriptFileName", required=True, \
		help='Name of file where script will be recorded')
	options = parser.parse_args();
	return options


def makeSummitCenteredFastaFromNarrowPeakScript(options):
	# Make a script that will predict the values for a new set of sequences
	narrowPeakFileNameListFile = open(options.narrowPeakFileNameListFileName)
	genomeFileNameListFile = open(options.genomeFileNameListFileName)
	chromSizesFileNameListFile = open(options.chromSizesFileNameListFileName)
	scriptFile = open(options.scriptFileName, 'w+')

	for narrowPeakFileNameStr, genomeFileNameStr, chromSizesFileNameStr in \
		zip(narrowPeakFileNameListFile, genomeFileNameListFile, \
			chromSizesFileNameListFile):
		# Iterate through the lines of the narrowPeak file list and the corresponding genomes and make a line in the script for each
		narrowPeakFileName = narrowPeakFileNameStr.strip()

		scriptFile.write(" ".join(["python", options.path + \
			"/makeSummitCenteredFastaFromNarrowPeak.py --narrowPeakFileName", \
			narrowPeakFileName, "--sequenceLength", \
			str(options.sequenceLength), "--genomeFileName", \
			genomeFileNameStr.strip(), "--chromSizesFileName", \
			chromSizesFileNameStr.strip()]))

		if options.createOptimalBed:
			# Add the createOptimalBed option
			scriptFile.write(" --createOptimalBed")
		if options.nameOnly:
                        # Add the createOptimalBed option
                        scriptFile.write(" --nameOnly")
		scriptFile.write("\n")

	narrowPeakFileNameListFile.close()
	genomeFileNameListFile.close()
	chromSizesFileNameListFile.close()
	scriptFile.close()


if __name__=="__main__":
	import sys
	options = parseArgument()
	makeSummitCenteredFastaFromNarrowPeakScript(options)
