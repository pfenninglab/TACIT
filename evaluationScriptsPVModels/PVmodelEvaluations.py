#usage: Evals_TESTset.py <model.h5>

import pygpu
import keras
import matplotlib
import sys
import os
import numpy as np
from keras.initializers import glorot_uniform
import theano
from Bio import SeqIO
import keras.backend as K
from keras.models import Sequential
from keras.layers import Dense, Conv2D, Flatten, Dropout, MaxPooling2D, Activation
from keras.regularizers import l2
from keras.optimizers import SGD
from keras.models import load_model
from keras import metrics
from sklearn.metrics import auc, roc_curve, precision_recall_curve, confusion_matrix
from sklearn.metrics import roc_auc_score
from keras.callbacks import Callback
from keras.callbacks import ModelCheckpoint
from sklearn.metrics import auc, roc_curve, precision_recall_curve, confusion_matrix
import sys

def onehot_seq(seq):
    letter_to_index =  {'A':0, 'a':0,
                        'C':1, 'c':1,
                        'G':2, 'g':2,
                        'T':3, 't':3}
    to_return = np.zeros((len(seq),4), dtype='int8')
    for idx,letter in enumerate(seq):
        if letter not in ['N','n']:
            to_return[idx,letter_to_index[letter]] = 1
    return to_return


def encode_sequence(fasta_pos, fasta_neg, shuffleOff = True):
    x_pos = np.array([onehot_seq(seq) for seq in SeqIO.parse(fasta_pos, "fasta") ] +
    [onehot_seq(seq.reverse_complement()) for seq in SeqIO.parse(fasta_pos, "fasta") ])
    x_neg = np.array([onehot_seq(seq) for seq in SeqIO.parse(fasta_neg, "fasta") ] +
    [onehot_seq(seq.reverse_complement()) for seq in SeqIO.parse(fasta_neg, "fasta") ])
    # concatenate positives and negatives
    x = np.expand_dims(np.concatenate((x_pos, x_neg)), axis=3)
    y = np.concatenate((np.ones(len(x_pos)),np.zeros(len(x_neg))))
    # need to shuffle order of training set for validation splitting last
    if not shuffleOff:
        indices = np.arange(y.shape[0])
        np.random.shuffle(indices)
        x = x[indices,:]
        y = y[indices]
    #
    return x, y

model = keras.models.load_model(sys.argv[1])

if os.path.exists("Eval1_mm10_TEST.fa") == False:
    sys.exit('ERROR: No file Eval1_mm10_TEST.fa')
if os.path.exists("Eval2_hg38_TEST.fa") == False:
    sys.exit('ERROR: No file Eval2_hg38_TEST.fa')
if os.path.exists("Eval3_hg38_TEST.fa") == False:
    sys.exit('ERROR: No file Eval3_hg38_TEST.fa')
if os.path.exists("Eval4_mm10_TEST.fa") == False:
    sys.exit('ERROR: No file Eval4_mm10_TEST.fa')
if os.path.exists("Eval5_mm10_TEST.fa") == False:
    sys.exit('ERROR: No file Eval5_mm10_TEST.fa')
if os.path.exists("Eval6_mm10_TEST.fa") == False:
    sys.exit('ERROR: No file Eval6_mm10_TEST.fa')
if os.path.exists("Eval7_mm10_TEST.fa") == False:
    sys.exit('ERROR: No file Eval7_mm10_TEST.fa')
if os.path.exists("Eval8_mm10_TEST.fa") == False:
    sys.exit('ERROR: No file Eval8_mm10_TEST.fa')
if os.path.exists("Eval9_mm10_TEST.fa") == False:
    sys.exit('ERROR: No file Eval9_mm10_TEST.fa')
if os.path.exists("Eval10_mm10_TEST.fa") == False:
    sys.exit('ERROR: No file Eval10_mm10_TEST.fa')



#multispecies evaluations (all restricted to TEST set for figure generation)

#Fig 1 Bar 1: + mouse&human positives, - mouse&human negatives
(x_eval12, y_eval12) = encode_sequence("/projects/pfenninggroup/mouseCxStr/NeuronSubtypeATAC/Zoonomia_CNN/multispecies_PV/FinalModelData/combined_pos_TEST.fa","/projects/pfenninggroup/mouseCxStr/NeuronSubtypeATAC/Zoonomia_CNN/multispecies_PV/FinalModelData/combined_neg_TEST.fa",shuffleOff=True)
y_pred = model.predict(x_eval12).ravel()
fpr_keras, tpr_keras, thresholds_keras = roc_curve(y_eval12, y_pred)
auc_12v = auc(fpr_keras, tpr_keras)
precision, recall, thresholds = precision_recall_curve(y_eval12, y_pred)
auprc_12v = auc(recall, precision)
y_predclass = np.rint(y_pred)
tn, fp, fn, tp = confusion_matrix(y_eval12, y_predclass).ravel()
acc_1v = tp/(tp+fn)
acc_2v = tn/(tn+fp)

#Fig 1 Bar 2: + mouse specific enhancers, - mouse non-enhancers whose human orthologs are enhancers
(x_eval14, y_eval14) = encode_sequence("Eval1_mm10_TEST.fa","Eval4_mm10_TEST.fa",shuffleOff=True)
y_pred = model.predict(x_eval14).ravel()
fpr_keras, tpr_keras, thresholds_keras = roc_curve(y_eval14, y_pred)
auc_14v = auc(fpr_keras, tpr_keras)
precision, recall, thresholds = precision_recall_curve(y_eval14, y_pred)
auprc_14v = auc(recall, precision)
y_predclass = np.rint(y_pred)
tn, fp, fn, tp = confusion_matrix(y_eval14, y_predclass).ravel()
acc_1v = tp/(tp+fn)
acc_1v = tn/(tn+fp)
# #pos > #neg, get aunpv
(x_eval41, y_eval41) = encode_sequence("Eval4_mm10_TEST.fa","Eval1_mm10_TEST.fa",shuffleOff=True)
y_pred41 = model.predict(x_eval41).ravel()
y_pred41 = y_pred41 * -1
precision, recall, thresholds = precision_recall_curve(y_eval41, y_pred41)
aunpv41 = auc(recall, precision)


#Fig 1 Bar 3: + human specific enchancers, - human non-enhancers whose mouse orthologs are enhancers
(x_eval32, y_eval32) = encode_sequence("Eval3_hg38_TEST.fa","Eval2_hg38_TEST.fa",shuffleOff=True)
y_pred = model.predict(x_eval32).ravel()
fpr_keras, tpr_keras, thresholds_keras = roc_curve(y_eval32, y_pred)
auc_32v = auc(fpr_keras, tpr_keras)
precision, recall, thresholds = precision_recall_curve(y_eval32, y_pred)
auprc_32v = auc(recall, precision)
y_predclass = np.rint(y_pred)
tn, fp, fn, tp = confusion_matrix(y_eval32, y_predclass).ravel()
acc_3v = tp/(tp+fn)
acc_2v = tn/(tn+fp)


#Fig 1 Bar 6: + shared liver/PV enhancers, - liver enhancers which are not PV enhancers
(x_eval910, y_eval910) = encode_sequence("Eval9_mm10_TEST.fa","Eval10_mm10_TEST.fa",shuffleOff=True)
y_pred = model.predict(x_eval910).ravel()
fpr_keras, tpr_keras, thresholds_keras = roc_curve(y_eval910, y_pred)
auc_910v = auc(fpr_keras, tpr_keras)
precision, recall, thresholds = precision_recall_curve(y_eval910, y_pred)
auprc_910v = auc(recall, precision)
y_predclass = np.rint(y_pred)
tn, fp, fn, tp = confusion_matrix(y_eval910, y_predclass).ravel()
acc_9v = tp/(tp+fn)
acc_10v = tn/(tn+fp)


#Fig 1 Bar 7: + shared PV & excitatory neuron enhancers, - excitatory neuron enhancers which are note PV enhancers
(x_eval67, y_eval67) = encode_sequence("Eval6_exc_MoHu_TEST.fa","Eval7exc_MoHu_TEST.fa",shuffleOff=True)
y_pred = model.predict(x_eval67).ravel()
fpr_keras, tpr_keras, thresholds_keras = roc_curve(y_eval67, y_pred)
auc_67v = auc(fpr_keras, tpr_keras)
precision, recall, thresholds = precision_recall_curve(y_eval67, y_pred)
auprc_67v = auc(recall, precision)
y_predclass = np.rint(y_pred)
tn, fp, fn, tp = confusion_matrix(y_eval67, y_predclass).ravel()
acc_6v = tp/(tp+fn)
acc_7v = tn/(tn+fp)


with open(str(sys.argv[1]) + 'evaluations1-10_TEST.txt', 'w') as f:
    f.write(str(auc_12v) + "\t" + str(auc_14v) + "\t" + str(auc_32v) + "\t" + str(auc_910v) + "\t" + str(auc_67v))
    f.write(str(auprc_12v) + "\t" + str(aunpv41) + "\t" + str(auprc_32v) + "\t" + str(auprc_910v) + "\t" + str(auprc_67v))
