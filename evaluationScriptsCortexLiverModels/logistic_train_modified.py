import sys
from math import log10
import numpy as np
from sklearn import linear_model as lm
from sklearn import metrics
from pickle import load, dump
from os.path import exists

#Arguments: Mode ("train" or not), fimo results, model file, predictions file
training = sys.argv[1].lower() == "train"
predicting = sys.argv[1].lower() == "predict"
preprocess = (exists(sys.argv[4]) and (predicting == False))
if training:
    motifs = set()
    lrModel = None
else:
    inModel = open(sys.argv[3], "rb")
    lrModel, sortedMotifs = load(inModel)
    motifs = set(sortedMotifs)
    inModel.close()
if preprocess:
    inData = open(sys.argv[4], "rb")
    peaks, features, labels, preds = load(inData)
    inData.close()
else:
    peaks = []
    if not predicting:
        # Evaluating the model, so define peak labels
        peakLabels = []
    peakDict = {}

            
    #motifs = set()
    fimoFile = open(sys.argv[2])
    fimoFile.readline()
    for line in fimoFile:
        #tokens = line.strip().split("\t")
        tokens = line.strip().split()
        if len(tokens) <= 1:
            continue
        motifName = tokens[1]
        peakName = tokens[2]
        if not predicting:
            # Evaluating the model, so define peak labels
            peakLabel = 1 if peakName.split("_")[-1] == "+" else 0
        pValue = float(tokens[-2])
        if peakName not in peakDict:
            peaks.append(peakName)
            if not predicting:
                # Evaluating the model, so define peak labels
                peakLabels.append(peakLabel)
            peakDict[peakName] = {}
        if motifName in peakDict[peakName]:
            peakDict[peakName][motifName] = min(pValue, peakDict[peakName][motifName])
        else:
            if motifName in motifs:
                peakDict[peakName][motifName] = pValue
            elif training:
                peakDict[peakName][motifName] = pValue
                motifs.add(motifName)    
                
    fimoFile.close()        
    print("Loaded data.")           

    if training:
        sortedMotifs = list(motifs)
        sortedMotifs.sort()
    features = np.zeros((len(peaks), len(sortedMotifs)))
    if not predicting:
        # Evaluating the model, so define labels
        labels = np.asarray(peakLabels)
    for i in range(len(peaks)):
        curPeakDict = peakDict[peaks[i]]
        for j in range(len(sortedMotifs)):
            pValue = curPeakDict.get(sortedMotifs[j])
            if pValue is not None:
                features[i][j] = -log10(pValue) - 3
    print("Processed data.")

if training:
    lrModel = lm.LogisticRegression(max_iter=1000, n_jobs=-1)
    lrModel.fit(features, labels)
    #print("Training accuracy", lrModel.score(features, labels))
    of = open(sys.argv[3], "wb")
    dump((lrModel, sortedMotifs), of, protocol=5)
    of.close()

if not predicting:
    # Evaluating the model instead of predicting on new sequences
    print("Accuracy", lrModel.score(features, labels))
if not preprocess:
    of = None
    if not predicting:
        # Have labels, so obtaining classes and writing predictions and labels to pkl file
        of = open(sys.argv[4], "wb")
        preds = lrModel.predict(features)
        dump((peaks, features, labels, preds), of, protocol=5)
    else:
        # Making predictions on new data, so obtaining probabilities and writing new peak name and prediction to a text file
        of = open(sys.argv[4], "w+")
        preds = lrModel.predict_proba(features)
        for peak, pred in zip(peaks, preds):
            of.write(peak + "\t" + str(pred[1]) + "\n")
    of.close()

if not predicting:
    # Evaluating the model instead of predicting on new sequences
    print("Performance of model:")
    print(metrics.classification_report(labels, preds, digits=4))
    print()
    print("ROC AUC:", metrics.roc_auc_score(labels, preds))
    p, r, t = metrics.precision_recall_curve(labels, preds)
    print("AUPRC:", metrics.auc(r, p))
    pNeg, rNeg, tNeg = metrics.precision_recall_curve(1 - labels, 1-preds)
    print("AUNPV-Spec.:", metrics.auc(rNeg, pNeg))
