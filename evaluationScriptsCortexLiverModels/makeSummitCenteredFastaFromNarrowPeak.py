import sys
import argparse
import os
import itertools
import numpy as np
import pybedtools as bt
from sequenceOperations import \
	createPositiveSetFromNarrowPeaks

"""
This script makes summit-centered fasta file for narrowPeak file and corresponding fasta and chromosome sizes file.
To run: python makeSummitCenteredFastaFromNarrowPeak.py [options]
"""

def parseArgument():
	# Parse the input
	parser = argparse.ArgumentParser(description = \
		"Make a summit-centered fasta from a narrowPeak file")
	parser.add_argument("--narrowPeakFileName", required=False, \
		help='Regions in narrowPeak format')
	parser.add_argument("--maxPeakLength", type=int, required=False, \
		default=None, help='Maximum number of bases per peak')
	parser.add_argument("--sequenceLength", type=int, required=False, \
		default=500, help='Number of bases in each sequence')
	parser.add_argument("--genomeFileName", required=False, \
		help='File with genome sequence')
	parser.add_argument("--chromSizesFileName", required=False, \
		help='File with list of chromosomes and their corresponding sizes')
	parser.add_argument("--chromEdgeDistLimit", type=int, required=False, default=0, \
                help='Distance from chromosome end for which peak will be considered, \
			need chromSizesFileName to use for 3 prime end of chormosomes')
	parser.add_argument("--chroms", action='append', required=False, \
		help='Chromosomes for which the file will be made')
	parser.add_argument("--createOptimalBed", action='store_true', required=False, \
                help='Remove peaks on unknown, alternative, \
			and random chromosomes before making predictions')
	parser.add_argument("--nameOnly", action='store_true', required=False, \
		help='Use the nameOnly option when creating the fasta file')
	options = parser.parse_args();
	return options

if __name__=="__main__":
	options = parseArgument()
	createPositiveSetFromNarrowPeaks(options.narrowPeakFileName, \
        	options.genomeFileName, \
                dataShape=(1,4,options.sequenceLength), \
                createOptimalBed=options.createOptimalBed, \
                createOptimalBedFilt=True, \
                maxPeakLength=options.maxPeakLength, \
                chroms=options.chroms, \
                chromSizesFileName=options.chromSizesFileName, \
                chromEdgeDistLimit=options.chromEdgeDistLimit, nameOnly=options.nameOnly)
	bt.helpers.cleanup()
