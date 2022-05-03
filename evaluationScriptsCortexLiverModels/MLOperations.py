import numpy as np
from sklearn.metrics import precision_recall_curve, recall_score, auc

def epsilonInsensitiveLoss(y_true, y_pred, epsilon=1.0):
	# Compute the epsilon-insensitive loss
	import keras.backend.theano_backend as K
	lossList = K.maximum(K.abs(y_true - y_pred) - epsilon, 0.)
	return K.mean(lossList)
	
def reweightedCrossentropyLoss(y_true, y_pred):
	# Compute the re-weighted cross-entropy loss, where everything with label -1 is weighted to 0
	import keras.backend.theano_backend as K
	import theano.tensor as T
	nonAmbig = (y_true > -0.5)
	cnts = nonAmbig.sum(axis=0, keepdims=True)*1.0
	assert (T.gt(cnts, 0.0))
	return K.mean(K.binary_crossentropy(y_pred*nonAmbig, y_true*nonAmbig)*y_true.shape[0]*1.0/cnts, axis=-1)
	
def taskweightedCrossentropyLoss(y_true, y_pred):
	# Compute the task-weighted cross-entropy loss, where every task is weighted by 1 - (fraction of non-ambiguous examples that are positive)
	# In addition, weight everything with label -1 to 0
	import keras.backend.theano_backend as K
	import theano.tensor as T
	nonAmbig = (y_true > -0.5)
	cnts = nonAmbig.sum(axis=0, keepdims=True)*1.0
	assert (T.gt(cnts, 0.0)) # Prevents division by 0
	pos = (y_true > 0.5)
	posCntsPerTask = pos.sum(axis=0, keepdims=True)*1.0
	assert (T.gt(posCntsPerTask, 0.0)) # Ensures that all tasks have at least 1 positive example
	weightsPerTask = 1.0 - (posCntsPerTask/cnts)
	weightsPerTaskRep = T.extra_ops.repeat(weightsPerTask, y_true.shape[0], axis=0)
	nonAmbigTimesWeightsPerTask = nonAmbig * weightsPerTaskRep
	normConst = nonAmbigTimesWeightsPerTask.sum(axis=0, keepdims=True).sum(axis=1)*1.0
	assert (T.gt(normConst, 0.0)) # Prevents division by 0
	return K.mean(K.binary_crossentropy(y_pred*nonAmbigTimesWeightsPerTask, y_true*nonAmbigTimesWeightsPerTask)*y_true.shape[0]*y_true.shape[1]*1.0/normConst, axis=-1)

# Compute the logistic loss for a single element
def elementLogisiticLoss(label, prediction):
	negLogisticLoss = (label*np.log(prediction)) + ((1 - label) * np.log(1 - prediction))
	return -negLogisticLoss

# Compute the specificity
def negativeAccuracy(labels, predictedLabels):
	labelsReverse = 1 - labels
	predictedLabelsReverse = 1 - predictedLabels
	return recall_score(labelsReverse, predictedLabelsReverse)
	
# Get the AUPRC
def auPRC(y_true, y_score):
	precision, recall, thresholds = precision_recall_curve(y_true, y_score)
	return auc(recall, precision)

# Get the recall at a specific FDR
def recallAtFdr(y_true, y_score, fdr_cutoff=0.05, returnCurve=False):
	precision, recall, thresholds = precision_recall_curve(y_true, y_score)
	precisionCutoff = 1 - fdr_cutoff
	cutoffIndex = 0
	for i in range(len(precision)):
		if precision[i] >= precisionCutoff:
			# At the precision cutoff
			cutoffIndex = i
			break
	if returnCurve:
		# Return the PR curve in addition to the recall at the desired FDR
		return [recall[cutoffIndex], precision, recall]
	return recall[cutoffIndex]
