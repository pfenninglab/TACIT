import csv
import os
import sys
import numpy as np
import pandas as pd
import pygpu
import keras
import matplotlib
from keras.models import load_model
from Bio import SeqIO

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

with open("tmpfFileList.txt", "r") as in1:
    fFileList = in1.read().splitlines()

with open("rows_peaks.txt", "r") as in2:
    rows_peaks = in2.read().splitlines()

model = keras.models.load_model("/projects/pfenninggroup/mouseCxStr/NeuronSubtypeATAC/Zoonomia_CNN/mouse_PV/models/modelPV3e.h5")


for f in fFileList:
    species = f[136:len(f)-3]
    print("scoring "+species)
    x_fwd = np.array([onehot_seq(seq) for seq in SeqIO.parse(f, "fasta") ] )
    xf = np.expand_dims((x_fwd), axis=3)
    fwd_pred = model.predict(xf).ravel()
    x_rev = np.array([onehot_seq(seq.reverse_complement()) for seq in SeqIO.parse(f, "fasta") ] )
    xr = np.expand_dims((x_rev), axis=3)
    rev_pred = model.predict(xr).ravel()
    avg_pred = np.mean([fwd_pred,rev_pred], axis=0)
    p = pd.DataFrame(index=rows_peaks,columns=[species])
    peakKeyFileName = "oldModel/"+species+"_peakKey.txt"
    with open(peakKeyFileName, "r") as in3:
        peakKey = in3.read().splitlines()
    for l in range(0,len(peakKey)):
        p.at[peakKey[l],species] = avg_pred[l]
    p.to_csv(r'pred_mouseReproduciblePeaks_'+species+'_modelPV3e.txt', sep="\t", mode='a', na_rep='NA')
