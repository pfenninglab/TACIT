import sys
import argparse
from Bio import SeqIO
from Bio.SeqIO import FastaIO

def parseArgument():
        # Parse the input
        parser=argparse.ArgumentParser(description=\
                        "Make the reverse complements of all sequences in a fasta file")
        parser.add_argument("--fastaFileName", required=True,\
                        help='Name of the original fasta file')
        parser.add_argument("--outputFileName", required=True,\
                        help='fasta file where the sequences that are the reverse complements will be written')
        options = parser.parse_args()
        return options

def reverseComplementFasta(options):
	# Make the reverse complements of all sequences in a fasta file
	seqRecordRCList = []
	outputFile = open(options.outputFileName, 'w+')
	fastaOut = FastaIO.FastaWriter(outputFile, wrap=None)
	for seqRecord in SeqIO.parse(options.fastaFileName,"fasta"):
		# Iterate through the sequences and get the reverse complement of each
		seqRecordRC = seqRecord.reverse_complement()
		seqRecordRC.id = seqRecord.id + "_RC"
		seqRecordRCList.append(seqRecordRC)
	fastaOut.write_file(seqRecordRCList)

if __name__ == "__main__":
        options = parseArgument()
        reverseComplementFasta(options)
