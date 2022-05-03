import sys
import argparse
import numpy as np
from scipy.stats import pearsonr, spearmanr
from keras.models import model_from_json
from sklearn.metrics import log_loss, roc_auc_score, accuracy_score, precision_recall_curve, auc, f1_score, recall_score, precision_score, \
	mean_squared_error, r2_score, roc_curve
from sklearn.linear_model import LogisticRegression as LR
from sequenceOperations import makeSequenceInputArrays, makeSequenceInputArraysNoLabels, \
	makeSequenceInputArraysFromDifferentialPeaks, convertFastaFileToSequencesFile, createPositiveSetFromNarrowPeaks
from MLOperations import negativeAccuracy, recallAtFdr, auPRC
#from auPRG_with_rpy2 import auPRG_R
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
	
def parseArgument():
	# Parse the input
	parser = argparse.ArgumentParser(description = "Predict the values for a new set of sequences and evalute the predictions")
	parser.add_argument("--architectureFileName", required=True, help='Model architecture file in json format')
	parser.add_argument("--weightsFileName", required=True, help='Model weights file in hdf5 format')
	parser.add_argument("--fastaFileName", required=False, help='Sequences file in fasta format')
	parser.add_argument("--negFastaFileName", required=False, help='Negative set sequences file in fasta format, assumes single-task model')
	parser.add_argument("--narrowPeakFileName", required=False, help='Regions in narrowPeak format')
	parser.add_argument("--negNarrowPeakFileName", required=False, help='Negative set regions in narrowPeak format, assumes single-task model')
	parser.add_argument("--labelsFileName", required=False, help='Labels/signals file in txt format')
	parser.add_argument("--label", type=int, action="append", required=False, default=None, help='Label that all of the sequences have for each task')
	parser.add_argument("--maxPeakLength", type=int, required=False, default=None, help='Maximum number of bases per peak')
	parser.add_argument("--sequenceLength", type=int, required=False, default=1000, help='Number of bases in each sequence')
	parser.add_argument("--perBaseTrackFileName", action='append', required=False, help='Per base track file name (like conservation, DNA shape, etc.')
	parser.add_argument("--predictedClassesFileName", default=None, required=False, help='File where predicted classes will be recorded')
	parser.add_argument("--predictedProbaFileName", default=None, required=False, help='File where predicted probabilities will be recorded')
	parser.add_argument("--differentialPeaks", action='store_true', required=False, help='Include if data is from output of DESeq2')
	parser.add_argument("--DESeq2OutputFileName", required=False, help='Output file from DESeq2 with data')
	parser.add_argument("--genomeFileName", required=False, help='File with genome sequence')
	parser.add_argument("--chromSizesFileName", required=False, help='File with chromosome sizes')
	parser.add_argument("--chromEdgeDistLimit", type=int, required=False, default=0, \
		help='Distance from chromosome end for which peak will be considered, need chromSizesFileName to use for 3 prime end of chormosomes')
	parser.add_argument("--backgroundFileName", required=False, help='File with the peaks or peak summits that went into DESeq2')
	parser.add_argument("--chroms", action='append', required=False, help='Chromosomes for which the prediction will happen')
	parser.add_argument("--maxFracNs", type=float, required=False, default=1.0, help='Maximum fraction of N\'s per sequence with prediction')
	parser.add_argument("--removeFastas", action='store_true', required=False, \
		help='Remove the fasta files that are created before making the numpy arrays')
	parser.add_argument("--swapLabels", action='store_true', required=False, \
		help='Swap the negative and positive labels for evaluation, only for differential peaks')
	parser.add_argument("--RC", action='store_true', required=False, help='Make sequence reverse complements')
	parser.add_argument("--logLabels", action='store_true', required=False, help='Include if labels (can be signals) should be log2ed')
	parser.add_argument("--multiMode", action='store_true', required=False, help='Include if model is multi-modal')
	parser.add_argument("--classification", action='store_true', required=False, help='Include if evaluating classification')
	parser.add_argument("--evaluateConsecutive", action='store_true', required=False, help='Include if binding and non-binding regions are consecutive')
	parser.add_argument("--recalibrateModel", action='store_true', required=False,\
		help='Recalibrate model before making predictions, implemented for classification only')
	parser.add_argument("--calibrationXTrFileName", required=False, help='npy file with calibration data')
	parser.add_argument("--calibrationYTrFileName", required=False, help='npy file with calibration labels')
	parser.add_argument("--plotRocPr", action='store_true', required=False, help='Plot the ROC and PR curves')
	parser.add_argument("--plotTitleSuffix", required=False, default="Validation-Set", help='Suffix of title for ROC and PR curve plots')
	parser.add_argument("--ROCFileName", required=False, help='Name of ROC plot file, should end in png')
	parser.add_argument("--PRFileName", required=False, help='Name of PR plot file, should end in png')
	parser.add_argument("--path", required=False, default="/srv/scratch/imk1/TFBindingPredictionProject/src/deeplearning/", \
		help='Path of deeplearning directory, should end with /')
	options = parser.parse_args();
	return options

def predictNewSequencesClassification(options, model, X, Y):
	print("Evaluating model...")
	predictedProba = model.predict_proba(X);
	predictedClasses = None;
	if Y.shape[1] == 1:
		predictedClasses = model.predict_classes(X);
	else:
		# There are multiple tasks, so the predict_classes method cannot be used (it takes the max. probability over the tasks)
                predictedClasses = np.int8(np.round(predictedProba));
	if options.recalibrateModel:
                # Re-calibrate the model with different data than what was used for training
                print("Calibrating the model")
                calibrator = LR()
		calibrationXTr = np.load(options.calibrationXTrFileName)
		calibrationYTr = np.load(options.calibrationYTrFileName)
                calibrationTrainPredictions = model.predict_proba(calibrationXTr)
                if len(calibrationTrainPredictions.shape) == 1:
                        # Expand the array
                        calibrationTrainPredictions = calibrationTrainPredictions[:,np.newaxis]
                calibrator.fit(calibrationTrainPredictions, calibrationYTr)
		for i in range(Y.shape[1]):
			# Iterate through the tasks and get the calibrated prediction for class 1 for each task
			predictedProba[:,i] = calibrator.predict_proba(predictedProba[:,i][:,np.newaxis])[:,1]
		predictedClasses = np.int8(np.round(predictedProba))
	print("Saving predictions (if files have been specified")
	if options.predictedClassesFileName != None:
		# Record the predicted classes
		np.savetxt(options.predictedClassesFileName, predictedClasses, fmt='%f')
	if options.predictedProbaFileName != None:
		# Record the predicted classes
		np.savetxt(options.predictedProbaFileName, predictedProba, fmt='%f')
	for i in range(Y.shape[1]):
		# Iterate through the tasks and get the results for each
		print("Results for task " + str(i))
		# accuracy
		acc = accuracy_score(Y[(Y[:,i] >= 0),i], predictedClasses[(Y[:,i] >= 0),i]);
		if options.label == None:
			# sensitivity
			sens = recall_score(Y[(Y[:,i] >= 0),i], predictedClasses[(Y[:,i] >= 0),i]);
			# specificity
			spec = negativeAccuracy(Y[(Y[:,i] >= 0),i], predictedClasses[(Y[:,i] >= 0),i]);
			# roc auc
			AUC = roc_auc_score(Y[(Y[:,i] >= 0),i], predictedProba[(Y[:,i] >= 0),i]);
			# precision
			precision = precision_score(Y[(Y[:,i] >= 0),i], predictedClasses[(Y[:,i] >= 0),i])
			# negative preictive value
			NPV = precision_score(1 - Y[(Y[:,i] >= 0),i], 1 - predictedClasses[(Y[:,i] >= 0),i])
			yTrue = Y[(Y[(Y[:,i] >= 0),i] >= 0)]
			yScore = predictedProba[(Y[(Y[:,i] >= 0),i] >= 0)]
			if np.sum(yTrue) > 0.5*yTrue.shape[0]:
				# Swap the negative and positive labels
				print("Majority of examples are positives, so evaluating for negative set")
				yTrue = 1 - yTrue
				yScore = 1 - yScore
			# auPRC
			AUPRC = auPRC(yTrue, yScore);
			# auNPV-Specificity
			AUNPVSpec = auPRC(1-yTrue, 1-yScore);
			# auPRG
			#AUPRG = auPRG_R(yTrue, yScore);
			# F1
			F1 = f1_score(Y[(Y[:,i] >= 0),i], predictedClasses[(Y[:,i] >= 0),i]);
			# Recall at FDR 0.05
			recallAt05FDR = recallAtFdr(yTrue, yScore, fdr_cutoff=0.05);
			# Recall at FDR 0.1
			[recallAt1FDR, precisionCurve, recallCurve] = recallAtFdr(yTrue, yScore, fdr_cutoff=0.1, returnCurve=True)
			# cross entropy loss
			loss = log_loss(Y[(Y[:,i] >= 0),i], predictedProba[(Y[:,i] >= 0),i]);
			# Print results
			#print("loss: %f, accuracy: %f, sensitivity: %f, specificity: %f, auc: %f, precision: %f, NPV: %f, auprc: %f, aunpvspec: %f, auprg: %f, F1: %f, recallAt0.05FDR: %f, recallAt0.1FDR: %f" \
			#	% (loss, acc, sens, spec, AUC, precision, NPV, AUPRC, AUNPVSpec, AUPRG, F1, recallAt05FDR, recallAt1FDR));
			print("loss: %f, accuracy: %f, sensitivity: %f, specificity: %f, auc: %f, precision: %f, NPV: %f, auprc: %f, aunpvspec: %f, F1: %f, recallAt0.05FDR: %f, recallAt0.1FDR: %f" \
                                % (loss, acc, sens, spec, AUC, precision, NPV, AUPRC, AUNPVSpec, F1, recallAt05FDR, recallAt1FDR));
		if options.plotRocPr:
			# Plot the ROC and PR curves
			assert(Y.shape[1] == 1)
			FPR, TPR, _ =  roc_curve(Y, predictedProba)
			lw = 6
			plt.plot(FPR, TPR, color='darkorange', lw=lw, label='ROC curve (area = %0.2f)' % AUC)
			plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
			plt.xlim([-0.05, 1.05])
			plt.ylim([-0.05, 1.05])
			plt.xlabel('False Positive Rate', fontsize=24)
			plt.ylabel('True Positive Rate', fontsize=24)
			plt.title('Receiver Operating Characteristic for ' + options.plotTitleSuffix, fontsize=24)
			plt.legend(loc="lower right", fontsize=20)
			plt.savefig(options.ROCFileName)
			plt.clf()
			plt.close()
			plt.plot(recallCurve, precisionCurve, color='darkorange', lw=lw, label='PR curve (area = %0.2f)' % AUPRC)
			plt.xlim([-0.05, 1.05])
			plt.ylim([-0.05, 1.05])
			plt.tick_params(labelsize=20)
			plt.xlabel('Recall')
			plt.ylabel('Precision')
			plt.title('Precision-Recall Curve for ' + options.plotTitleSuffix)
			plt.legend(loc="lower right", fontsize=20)
			plt.savefig(options.PRFileName)
			plt.clf()
			plt.close()
		else:
			# Print only the loss and the accuracy
			print("accuracy: %f" % (acc));

def predictNewSequencesRegression(model, X, Y):
	print("Evaluating model...")
	predictions = model.predict(X);
	if len(Y.shape) == 1:
		# Expand Y into a 2D array
		Y = Y[:,np.newaxis]
	for i in range(Y.shape[1]):
		# Iterate through the tasks and get the results for each task
		print ("Results for task " + str(i))
		# r^2
		rSquared = r2_score(Y[:,i], predictions[:,i])
		# Spearman correlation
		spearCorr = spearmanr(Y[:,i], predictions[:,i])
		# MSE loss
		loss = mean_squared_error(Y[:,i], predictions[:,i])
		# Pearson correlation
		pCorr = pearsonr(Y[:,i], predictions[:,i])
		# Print results
		print("loss: %f, r^2: %f, PearsonCorrelation: %f, PearsonCorrelationPVal: %f, SpearmanCorrelation: %f, SpearmanCorrelationPVal: %f" % \
			(loss, rSquared, pCorr[0], pCorr[1], spearCorr[0], spearCorr[1]));
		
def predictNewSequences(options):
	# Predict the values for a new set of sequences and evalute the predictions
	model = model_from_json(open(options.architectureFileName).read())
	print("Loading model")
	model.load_weights(options.weightsFileName)
	print("Model has been loaded.")
	X = np.array(1)
	Y = np.array(1)
	sequencesFileName = None
	numSequences = None
	negSequencesFileName = None
	numNegSequences = None
	if options.fastaFileName != None:
		# The sequences were given in fasta format, so convert their format
		sequencesFileName, numSequences, sequenceIDs = convertFastaFileToSequencesFile(options.fastaFileName)
		if options.negFastaFileName != None:
			# The fasta file with the sequences in the negative set have been provided in a separate file, so get the sequences of the negative set
			negSequencesFileName, numNegSequences, negSequenceIDs = convertFastaFileToSequencesFile(options.negFastaFileName)
	elif options.narrowPeakFileName != None:
		# Only the peaks have been provided, so get the sequences from the peak summits +/-
		_, _, _, _, _, positiveFastaFileName =\
			createPositiveSetFromNarrowPeaks(options.narrowPeakFileName, options.genomeFileName, \
				dataShape=(1,4,options.sequenceLength), createOptimalBed=False, createOptimalBedFilt=True, \
				maxPeakLength=options.maxPeakLength, chroms=options.chroms, chromSizesFileName=options.chromSizesFileName, \
				chromEdgeDistLimit=options.chromEdgeDistLimit)
		sequencesFileName, numSequences, sequenceIDs = convertFastaFileToSequencesFile(positiveFastaFileName)
		print ("The number of sequences is:  " + str(numSequences))
		if options.negNarrowPeakFileName != None:
			# The peaks from the negative set have been provided in a sepearte file, so get the sequences of the negative set
			_, _, _, _, _, negativeFastaFileName =\
                        	createPositiveSetFromNarrowPeaks(options.negNarrowPeakFileName, options.genomeFileName, \
                                	dataShape=(1,4,options.sequenceLength), createOptimalBed=False, createOptimalBedFilt=True, \
                                	maxPeakLength=options.maxPeakLength, chroms=options.chroms, chromSizesFileName=options.chromSizesFileName, \
					chromEdgeDistLimit=options.chromEdgeDistLimit)
			negSequencesFileName, numNegSequences, negSequenceIDs = convertFastaFileToSequencesFile(negativeFastaFileName)
			print ("The number of sequences in the negative set is:  " + str(numNegSequences))
	if options.differentialPeaks:
		assert (maxFracNs == 1.0)
		X, Y, _, _, _, _ = makeSequenceInputArraysFromDifferentialPeaks(options.DESeq2OutputFileName, options.genomeFileName, \
			options.backgroundFileName, (1,4,options.sequenceLength), createOptimalBed=False, backgroundSummitPresent=False, \
			backgroundSummitOnly=True, createModelDir=False, chroms=options.chroms, bigWigFileNames=[], \
			multiMode=options.multiMode, streamData=False, dataFileName="", RC=options.RC, removeFastas=options.removeFastas, \
			strictNegativeSet=False, fcCutoff=1, swapLabels=options.swapLabels, chromEdgeDistLimit=options.chromEdgeDistLimit)
	else:
		# The sequence is haploid
		if (options.label == None) and ((options.negNarrowPeakFileName == None) and (options.negFastaFileName == None)):
			# The labels are in the file of labels
			X, Y =\
				makeSequenceInputArrays(sequencesFileName, options.labelsFileName, (1,4,options.sequenceLength), logLabels=options.logLabels, \
					perBaseTrackFileNames=options.perBaseTrackFileName, multiMode=options.multiMode, RC=options.RC, \
					maxFracNs=options.maxFracNs)
		elif options.label != None:
			# All sequences share the same label
			assert (options.negNarrowPeakFileName == None)
			X, _ =\
				makeSequenceInputArraysNoLabels(sequencesFileName, (1,4,options.sequenceLength), numSequences, \
					perBaseTrackFileNames=options.perBaseTrackFileName, multiMode=options.multiMode, maxFracNs=options.maxFracNs)
			Y = np.zeros((X.shape[0], len(options.label)))
			for i in range(len(options.label)):
				# Iterate through the tasks and set the label for each task
				Y[:,i] = options.label[i]
		else:
			assert ((options.negNarrowPeakFileName != None) or (options.negFastaFileName != None))
			assert (options.classification)
			XPos, _ =\
                                makeSequenceInputArraysNoLabels(sequencesFileName, (1,4,options.sequenceLength), numSequences, \
                                        perBaseTrackFileNames=options.perBaseTrackFileName, multiMode=options.multiMode, maxFracNs=options.maxFracNs)
			XNeg, _ =\
                                makeSequenceInputArraysNoLabels(negSequencesFileName, (1,4,options.sequenceLength), numNegSequences, \
                                        perBaseTrackFileNames=options.perBaseTrackFileName, multiMode=options.multiMode, maxFracNs=options.maxFracNs)
			YPos = np.ones((XPos.shape[0], 1))
			YNeg = np.zeros((XNeg.shape[0], 1))
			X = np.vstack((XPos, XNeg))
			Y = np.vstack((YPos, YNeg))
	if options.classification:
		# Evaluate the classification model
		predictNewSequencesClassification(options, model, X, Y)
	else:
		# Evaluate the regression model
		predictNewSequencesRegression(model, X, Y)

if __name__=="__main__":
	options = parseArgument()
	predictNewSequences(options)
