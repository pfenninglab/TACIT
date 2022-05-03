import sys
import argparse
import gzip
import os

def parseArgument():
        # Parse the input
        parser=argparse.ArgumentParser(description=\
                        "Create a file with the predictions for each peak across all of the species")
        parser.add_argument("--predictionsFileNameListFileName", required=True,\
                        help='List of predictions files with the same peaks from different species, first column is peak names, 2nd column is predictions')
	parser.add_argument("--baseSpeciesLine", type=int, required=False, default=0,\
                        help='Number of line with the species from which open chromatin data was obtained, zero-indexed')
        parser.add_argument("--outputFileName", required=True, \
                        help='Name of file where the predictions will be recorded')
        options = parser.parse_args()
        return options

def gatherPeakPredictionsAcrossSpecies(options):
	# Create a file with the predictions for each peak across all of the species
	predictionsFileNameListFile = open(options.predictionsFileNameListFileName)
	peakDictList = []
	outputFile = open(options.outputFileName, 'w+')
	for predictionsFileNameLine in predictionsFileNameListFile:
		# Iterate through the bed files and create a dictionary and output file for each
		predictionsFileName = predictionsFileNameLine.strip()
		peakDict = {}
		if os.path.exists(predictionsFileName):
			# There are peak orthologs in the species corresponding to the current file
			predictionsFile = open(predictionsFileName)
			for line in predictionsFile:
				# Iterate through the lines of the bed file and add an entry to the dictionary for each
				lineElements = line.strip().split("\t")
				if len(lineElements) < 2:
					print (predictionsFileName + " is in the wrong format")
				assert(len(lineElements) >= 2)
				peakName = lineElements[0]
				# If the prediction for the peak has already been recorded (current prediction is for the RC), store the mean
				peakDict[peakName] =\
					((peakDict[peakName] + float(lineElements[1]))/2.0) if peakName in peakDict else \
						float(lineElements[1])
			predictionsFile.close()
		else:
			# There are no peak orthologs in the species corresponding to the current file
			print(predictionsFileName + " does not exist.")
		peakDictList.append(peakDict)
	predictionsFileNameListFile.close()
	for peak in peakDictList[options.baseSpeciesLine].keys():
		# Iterate through the peaks from the first bed file and get the predictions for each
		outputFile.write(peak)
		for peakDict in peakDictList:
			# Iterate through the species and the peak predictions if there is an ortholog in the current species, otherwise record -1
			outputFile.write("\t" + str(peakDict[peak]) \
				if peak in peakDict.keys() else "\t" + "-1")
		outputFile.write("\n")
	outputFile.close()

if __name__=="__main__":
        options = parseArgument()
        gatherPeakPredictionsAcrossSpecies(options)
