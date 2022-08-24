from tensorflow.keras.metrics import AUC, TruePositives, FalsePositives, TrueNegatives, FalseNegatives, Precision, Recall
from tensorflow.keras.models import load_model
import tensorflow as tf
import tensorflow.keras.backend as K
import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import argparse, sys, os
from tqdm import tqdm
from Bio import SeqIO
print(tf.__version__)
print(tf.keras.__version__)

# set seed
from numpy.random import seed
seed(1)
from tensorflow.random import set_seed
set_seed(2)

from sklearn.metrics import *
import matplotlib.pyplot as plt

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

def macro_f1(y, y_hat, thresh=0.5):
  y_pred = tf.cast(tf.greater(y_hat, thresh), tf.float32)
  tp = tf.cast(tf.math.count_nonzero(y_pred * y, axis=0), tf.float32)
  fp = tf.cast(tf.math.count_nonzero(y_pred * (1 - y), axis=0), tf.float32)
  fn = tf.cast(tf.math.count_nonzero((1 - y_pred) * y, axis=0), tf.float32)
  f1 = 2*tp / (2*tp + fn + fp + 1e-16)
  macro_f1 = tf.reduce_mean(f1, axis=-1)
  return macro_f1

def get_input(path):
    x = list()
    for seq in tqdm(SeqIO.parse(path, "fasta")):
      x.append(onehot_seq(seq))
      x.append(onehot_seq(seq.reverse_complement()))
    return np.expand_dims(np.array(x), axis=3)

def calc_metrics(model, x, y):
  # num pos, num neg, auroc, auprc, majority class auprc/aunpv spec, fraction pos/neg, sensitivity, specificity, precision, NPV
  num_pos = int(np.sum(y == 1)/2)
  num_neg = int(np.sum(y == 0)/2)
  y_pred = model.predict(x).ravel()

  precision_all, recall_all, thresholds = precision_recall_curve(y, y_pred)
  auPRC = auc(recall_all, precision_all)
  fraction = min(len(np.argwhere(y == 1)), len(np.argwhere(y == 0)))/(len(y))
  rounded_preds = y_pred.round()
  (tn, fp, fn, tp) = confusion_matrix(y, rounded_preds).ravel()
  sensitivity = tp/(tp+fn)
  specificity = tn/(tn+fp)

  precision = tp/(tp+fp)
  npv = tn/(tn+fn)
  swap_labels = np.where(y == 1, 0, 1)
  swap_rounded_preds = np.where(rounded_preds == 1, 0, 1)
  swap_precision_all, swap_recall_all, swap_thresholds = precision_recall_curve(swap_labels, swap_rounded_preds)
  aunpv_spec = auc(swap_recall_all, swap_precision_all)
  maj_prc_npv = aunpv_spec if num_neg >= num_pos else auPRC
  uPRC = auPRC if num_neg >= num_pos else aunpv_spec
  return num_pos, num_neg, roc_auc_score(y, y_pred), auPRC, maj_prc_npv, fraction, sensitivity, specificity, precision, npv

def print_out(model, xtest, ytest):
  num_pos, num_neg, auroc, auprc, maj_prc_npv, fraction, sensitivity, specificity, precision, npv = calc_metrics(model, xtest, ytest)
  print("positives", num_pos)
  print("negatives", num_neg)
  print("auroc", auroc)
  print("auprc", auprc)
  print("majority auprc/aunpv-spec", maj_prc_npv)
  print("fraction", fraction)
  print("sensitivity", sensitivity)
  print("specificity", specificity)
  print("precision", precision)
  print("npv", npv)

def main_fig(model):
  print("loading data...")
  mouse_data="/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/mouse/GSE146897_WTMouse_ret_ATAC/"
  human_data="/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/human/GSE137311/"
  eval1 = mouse_data+"mouse_specific_mm10.fa"
  eval2 = mouse_data+"mouse_specific_hg38.fa"
  eval3 = human_data+"human_specific_hg38.fa"
  eval4 = human_data+"human_specific_mm10.fa"
  eval5 = mouse_data+"mm10_ret_nooverlap_liver.fa"
  eval6 = mouse_data+"mm10_ret_overlap_liver.fa"
  eval7 = mouse_data+"mm10_ret_noTSS_filtered_500bp_neg_overlap_liver.fa"

  xtestpos = get_input("/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/mouse_human4/pos_TEST.fa")
  ytestpos = np.ones(len(xtestpos))
  xtestneg = get_input("/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/mouse_human4/neg_TEST.fa")
  ytestneg = np.zeros(len(xtestneg))
  xtest = np.concatenate((xtestpos, xtestneg))
  ytest = np.concatenate((ytestpos, ytestneg))

  x_eval1 = get_input(eval1)
  y_eval1 = np.ones(len(x_eval1))
  x_eval4 = get_input(eval4)
  y_eval4 = np.zeros(len(x_eval4))
  x_eval_mouse = np.concatenate((x_eval1, x_eval4))
  y_eval_mouse = np.concatenate((y_eval1, y_eval4))
  x_eval2 = get_input(eval2)
  y_eval2 = np.zeros(len(x_eval2))
  x_eval3= get_input(eval3)
  y_eval3 = np.ones(len(x_eval3))
  x_eval_human = np.concatenate((x_eval3, x_eval2))
  y_eval_human = np.concatenate((y_eval3, y_eval2))


  x_eval5 = get_input(eval5)
  y_eval5 = np.ones(len(x_eval5))
  x_eval6 = get_input(eval6)
  y_eval6 = np.ones(len(x_eval6))
  x_eval7 = get_input(eval7)
  y_eval7 = np.zeros(len(x_eval7))
  x_eval_tissue = np.concatenate((x_eval5, x_eval6, x_eval7))
  y_eval_tissue = np.concatenate((y_eval5, y_eval6, y_eval7))

  print("bar1")
  print_out(model, xtest, ytest)

  print("bar2")
  print_out(model, x_eval_mouse, y_eval_mouse)

  print("bar3")
  print_out(model, x_eval_human, y_eval_human)

  print("bar6")
  print_out(model, x_eval_tissue, y_eval_tissue)

def supp_fig1(model):

  print("loading data...")
  mouse_data="/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/mouse/GSE146897_WTMouse_ret_ATAC/"
  human_data="/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/human/GSE137311/"
  eval1 = mouse_data+"mouse_specific_mm10.fa"
  eval2 = mouse_data+"mouse_specific_hg38.fa"
  eval3 = human_data+"human_specific_hg38.fa"
  eval4 = human_data+"human_specific_mm10.fa"

  x_eval1 = get_input(eval1)
  y_eval1 = np.ones(len(x_eval1))
  x_eval4 = get_input(eval4)
  y_eval4 = np.zeros(len(x_eval4))
  x_eval_mouse = np.concatenate((x_eval1, x_eval4))
  y_eval_mouse = np.concatenate((y_eval1, y_eval4))
  x_eval2 = get_input(eval2)
  y_eval2 = np.zeros(len(x_eval2))
  x_eval3= get_input(eval3)
  y_eval3 = np.ones(len(x_eval3))
  x_eval_human = np.concatenate((x_eval3, x_eval2))
  y_eval_human = np.concatenate((y_eval3, y_eval2))

  xmousetestpos = get_input("/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/mouse4_v3/pos_TEST.fa")
  ymousetestpos = np.ones(len(xmousetestpos))
  xmousetestneg = get_input("/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/mouse4_v3/neg_TEST.fa")
  ymousetestneg = np.zeros(len(xmousetestneg))
  xmousetest = np.concatenate((xmousetestpos, xmousetestneg))
  ymousetest = np.concatenate((ymousetestpos, ymousetestneg))

  xtestpos = get_input("/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/mouse_human4/pos_TEST.fa")
  ytestpos = np.ones(len(xtestpos))
  xtestneg = get_input("/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/mouse_human4/neg_TEST.fa")
  ytestneg = np.zeros(len(xtestneg))
  xtest = np.concatenate((xtestpos, xtestneg))
  ytest = np.concatenate((ytestpos, ytestneg))

  print("bar1")
  #print(len(xmousetestpos), len(xmousetestneg))
  print_out(model, xmousetest, ymousetest)

  print("bar2")
  #print(len(xtestpos), len(xtestneg))
  print_out(model, xtest, ytest)

  print("bar3")
  #print(len(x_eval1), len(x_eval4))
  print_out(model, x_eval_mouse, y_eval_mouse)

  print("bar4")
  #print(len(x_eval2), len(x_eval3))
  print_out(model, x_eval_human, y_eval_human)


def supp_fig2(model):
  print("loading data...")
  mouse_data="/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/mouse/GSE146897_WTMouse_ret_ATAC/"
  human_data="/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/human/GSE137311/"
  eval1 = mouse_data+"mouse_specific_mm10.fa"
  eval2 = mouse_data+"mouse_specific_hg38.fa"
  eval3 = human_data+"human_specific_hg38.fa"
  eval4 = human_data+"human_specific_mm10.fa"
  eval5 = mouse_data+"mm10_ret_nooverlap_liver.fa"
  eval6 = mouse_data+"mm10_ret_overlap_liver.fa"
  eval7 = mouse_data+"mm10_ret_noTSS_filtered_500bp_neg_overlap_liver.fa"
  x_eval5 = get_input(eval5)
  y_eval5 = np.ones(len(x_eval5))
  x_eval6 = get_input(eval6)
  y_eval6 = np.ones(len(x_eval6))
  x_eval7 = get_input(eval7)
  y_eval7 = np.zeros(len(x_eval7))
  x_eval_tissue = np.concatenate((x_eval6, x_eval7))
  y_eval_tissue = np.concatenate((y_eval6, y_eval7))

  print("bar1")
  #print(np.sum(y_eval_tissue == 1), np.sum(y_eval_tissue == 0))
  print_out(model, x_eval_tissue, y_eval_tissue)
  #sys.exit()

  eval5 = mouse_data+"mm10_ret_nooverlap_cortex.fa"
  eval6 = mouse_data+"mm10_ret_overlap_cortex.fa"
  eval7 = mouse_data+"mm10_ret_noTSS_filtered_500bp_neg_overlap_cortex.fa"
  x_eval5 = get_input(eval5)
  y_eval5 = np.ones(len(x_eval5))
  x_eval6 = get_input(eval6)
  y_eval6 = np.ones(len(x_eval6))
  x_eval7 = get_input(eval7)
  y_eval7 = np.zeros(len(x_eval7))
  x_eval_tissue = np.concatenate((x_eval6, x_eval7))
  y_eval_tissue = np.concatenate((y_eval6, y_eval7))

  print("bar2")
  print_out(model, x_eval_tissue, y_eval_tissue)


def species_specific_metrics(model):
  humanposx = get_input("/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/human4/pos_TEST.fa")
  humanposy= np.ones(len(humanposx))
  humannegx = get_input("/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/human4/neg_TEST.fa")
  humannegy= np.zeros(len(humannegx))
  humantestx = np.concatenate((humanposx, humannegx))
  humantesty = np.concatenate((humanposy, humannegy))
  mouseposx = get_input("/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/mouse4_v3/pos_TEST.fa")
  mouseposy = np.ones((len(mouseposx)))
  mousenegx = get_input("/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/mouse4_v3/neg_TEST.fa")
  mousenegy = np.zeros((len(mousenegx)))
  mousetestx = np.concatenate((mouseposx, mousenegx))
  mousetesty = np.concatenate((mouseposy, mousenegy))
  #print("mouse metrics")
  print_out(model, mousetestx, mousetesty)
  print_out(model, humantestx, humantesty)

def new_project(model):
  mouse_data="/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/mouse/GSE146897_WTMouse_ret_ATAC/"
  human_data="/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/human/GSE137311/"
  eval2 = mouse_data+"mouse_specific_hg38.fa"
  eval3 = human_data+"human_specific_hg38.fa"
  x_eval2 = get_input(eval2)
  y_eval2 = np.zeros(len(x_eval2))
  x_eval3= get_input(eval3)
  y_eval3 = np.ones(len(x_eval3))
  y_pred = model.predict(x_eval2).ravel()
  np.save("/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/mouse4_v3/preds/human_specific_nonenh.npy", y_pred, allow_pickle=True)
  y_pred = model.predict(x_eval3).ravel()
  np.save("/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/mouse4_v3/preds/human_specific_enh.npy", y_pred, allow_pickle=True)


def train_validation_performance(model):
  postrain = "/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/mouse4_v3/pos_TRAINING.fa"
  negtrain = "/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/mouse4_v3/neg_TRAINING.fa"
  posvalid = "/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/mouse4_v3/pos_VALIDATION.fa"
  negvalid = "/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/mouse4_v3/neg_VALIDATION.fa"
  postrain = get_input(postrain)
  negtrain = get_input(negtrain)
  posvalid = get_input(posvalid)
  negvalid = get_input(negvalid)
  xtrain = np.concatenate((postrain, negtrain))
  ytrain = np.concatenate((np.ones(len(postrain)), np.zeros(len(negtrain))))
  print("training")
  print_out(model, xtrain, ytrain)
  xvalid = np.concatenate((posvalid, negvalid))
  yvalid = np.concatenate((np.ones(len(posvalid)), np.zeros(len(negvalid))))
  print("validation")
  print_out(model, xvalid, yvalid)



def main():
  print("loading model...")
  #multispeciesmodel = load_model("/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/mouse_human4/best_evals.h5", custom_objects={"macro_f1": macro_f1})
  mouseonly_model = load_model("/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/mouse4_heinit2.h5", custom_objects={"macro_f1": macro_f1})

  #main_fig(multispeciesmodel)
  #species_specific_metrics(multispeciesmodel)
  #supp_fig1(mouseonly_model)
  #supp_fig2(mouseonly_model)

  #new_project(mouseonly_model)
  train_validation_performance(mouseonly_model)



if __name__ == "__main__":
  main()
