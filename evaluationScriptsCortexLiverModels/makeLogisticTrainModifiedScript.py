import sys
import argparse

def parseArgument():
        # Parse the input
        parser = argparse.ArgumentParser(description =\
                "Make a script that will make logistic regression model predictions for a list of FIMO feature files")
        parser.add_argument("--FIMOFileNameListFileName", \
                required=True, \
                help='File with list of FIMO feature output file names')
        parser.add_argument("--modelFileName", \
                help='Name of file with logistic regression model in pkl format')
        parser.add_argument("--outputFileNameSuffix", required=False, default="_LRPredictions.txt", \
                help='Suffix of output file name for predictions file')
        parser.add_argument("--path", required=True, \
                help='Path to logistic_train_modified.py')
        parser.add_argument("--scriptFileName", required=True, \
                help='Name of file where script will be recorded')
        options = parser.parse_args();
        return options

def makeLogisticTrainModifiedScript(options):
	# Make a script that will make logistic regression model predictions for a list of FIMO feature files
	FIMOFileNameListFile = open(options.FIMOFileNameListFileName)
	scriptFile = open(options.scriptFileName, 'w+')
	fullFunctionName = options.path + "/logistic_train_modified.py"
	for line in FIMOFileNameListFile:
		# Iterate through the FIMO output file names and make a line in the script for each
		FIMOFileName = line.strip()
		FIMOFileNameElements = FIMOFileName.split(".")
		outputFileName = ".".join(FIMOFileNameElements[0:-1]) + options.outputFileNameSuffix
		scriptFile.write(" ".join(["python3", fullFunctionName, "predict", FIMOFileName, options.modelFileName, outputFileName]) + "\n")
	FIMOFileNameListFile.close()
	scriptFile.close()

if __name__=="__main__":
        import sys
        options = parseArgument()
        makeLogisticTrainModifiedScript(options)
