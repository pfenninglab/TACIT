import sys
import argparse
from itertools import izip
import numpy as np
from sklearn.metrics import log_loss, roc_auc_score, accuracy_score, f1_score, recall_score, precision_score, roc_curve
from MLOperations import negativeAccuracy, recallAtFdr, auPRC
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def parseArgument():
        # Parse the input
        parser = argparse.ArgumentParser(description = "Evaluate a set of predictions")
	parser.add_argument("--posPredFileName", action='append', required=True, help='File of predictions for positives')
	parser.add_argument("--posHeaderNum", action="append", type=int, required=True, help='Number of headers for each positive predictions file')
	parser.add_argument("--posPredCol", action='append', type=int, required=True, help='Column number of predictions for each positive predictions file')
	parser.add_argument("--negPredFileName", action='append', required=True, help='File of predictions for negatives')
        parser.add_argument("--negHeaderNum", action='append', type=int, required=True, help='Number of headers for each negative predictions file')
        parser.add_argument("--negPredCol", action='append', type=int, required=True, help='Column number of predictions for each negative predictions file')
        parser.add_argument("--plotRocPr", action='store_true', required=False, help='Plot the ROC and PR curves')
	parser.add_argument("--plotTitleSuffix", required=False, default="Test-Set", help='Suffix of title for ROC and PR curve plots')
        parser.add_argument("--ROCFileName", required=False, help='Name of ROC plot file, should end in png')
        parser.add_argument("--PRFileName", required=False, help='Name of PR plot file, should end in png')
        options = parser.parse_args();
        return options

def evaluateSequencePredictions(options):
	# Evaluate a set of predictions
	posPredictions = np.empty([0])
	for ppfn, ppn, ppc in izip(options.posPredFileName, options.posHeaderNum, options.posPredCol):
		# Iterate through the positive predictions and add each to the array
		posPred = np.loadtxt(ppfn, skiprows=ppn, usecols=ppc)
		posPredictions = np.hstack((posPredictions, posPred))
	negPredictions = np.empty([0])
        for npfn, npn, npc in izip(options.negPredFileName, options.negHeaderNum, options.negPredCol):
                # Iterate through the positive predictions and add each to the array
                negPred = np.loadtxt(npfn, skiprows=npn, usecols=npc)
                negPredictions = np.hstack((negPredictions, negPred))
	posLabels = np.ones((posPredictions.shape[0],), dtype=int)
	negLabels = np.zeros((negPredictions.shape[0],), dtype=int)
	predictions = np.hstack((posPredictions, negPredictions))
	predictedClasses = predictions.round()
	labels = np.hstack((posLabels, negLabels))
	print("Evaluating predictions")
	# accuracy
        acc = accuracy_score(labels, predictedClasses);
	# sensitivity
        sens = recall_score(labels, predictedClasses);
        # specificity
        spec = negativeAccuracy(labels, predictedClasses);
        # roc auc
        AUC = roc_auc_score(labels, predictions);
        # precision
        precision = precision_score(labels, predictedClasses)
        # negative preictive value
        NPV = precision_score(1 - labels, 1 - predictedClasses)
        yTrue = labels
        yScore = predictions
        if np.sum(labels) > 0.5*labels.shape[0]:
        	# Swap the negative and positive labels
                print("Majority of examples are positives, so evaluating for negative set")
                yTrue = 1 - yTrue
                yScore = 1 - yScore
        # auPRC
        AUPRC = auPRC(yTrue, yScore);
        # auNPV-Specificity
        AUNPVSpec = auPRC(1-yTrue, 1-yScore);
        # F1
        F1 = f1_score(labels, predictedClasses);
        # Recall at FDR 0.05
        recallAt05FDR = recallAtFdr(yTrue, yScore, fdr_cutoff=0.05);
        # Recall at FDR 0.1
        [recallAt1FDR, precisionCurve, recallCurve] = recallAtFdr(yTrue, yScore, fdr_cutoff=0.1, returnCurve=True)
        # cross entropy loss
        loss = log_loss(labels, predictions);
        # Print results
        print("loss: %f, accuracy: %f, sensitivity: %f, specificity: %f, auc: %f, precision: %f, NPV: %f, auprc: %f, aunpvspec: %f, F1: %f, recallAt0.05FDR: %f, recallAt0.1FDR: %f" \
        	% (loss, acc, sens, spec, AUC, precision, NPV, AUPRC, AUNPVSpec, F1, recallAt05FDR, recallAt1FDR));
	if options.plotRocPr:
        	# Plot the ROC and PR curves
                FPR, TPR, _ =  roc_curve(labels, predictions)
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

if __name__=="__main__":
        options = parseArgument()
        evaluateSequencePredictions(options)
