from tensorflow.keras.metrics import AUC, TruePositives, FalsePositives, TrueNegatives, FalseNegatives, Precision, Recall
from tensorflow.keras.models import load_model
import tensorflow as tf
import tensorflow.keras.backend as K
print(tf.__version__)
print(tf.keras.__version__)
from sklearn import linear_model
from sklearn.metrics import r2_score
from sklearn.metrics.pairwise import cosine_similarity
from scipy.stats import pearsonr, spearmanr
from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
import argparse, sys, os, glob
from tqdm import tqdm
from Bio import SeqIO
import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns

# set seed
from numpy.random import seed
seed(1)
from tensorflow.random import set_seed
set_seed(2)


def macro_f1(y, y_hat, thresh=0.5):
    """Compute the macro F1-score on a batch of observations (average F1 across labels)
    Args:
        y: (int32 Tensor)
            labels array of shape (BATCH_SIZE, N_LABELS)
        y_hat: (float32 Tensor)
            probability matrix from forward propagation of shape (BATCH_SIZE, N_LABELS)
        thresh: (float)
            probability value above which we predict positive

    Return
    ------
        macro_f1: (scalar Tensor)
            value of macro F1 for the batch
    """
    y_pred = tf.cast(tf.greater(y_hat, thresh), tf.float32)
    tp = tf.cast(tf.math.count_nonzero(y_pred * y, axis=0), tf.float32)
    fp = tf.cast(tf.math.count_nonzero(y_pred * (1 - y), axis=0), tf.float32)
    fn = tf.cast(tf.math.count_nonzero((1 - y_pred) * y, axis=0), tf.float32)
    f1 = 2*tp / (2*tp + fn + fp + 1e-16)
    macro_f1 = tf.reduce_mean(f1, axis=-1)
    return macro_f1


def onehot_seq(seq):
    '''one-hot-encode a DNA sequence

    Input
    -----
    seq: (string)
        DNA sequence

    Return
    ------
    to_return: (np.array, dim=(len(seq), 4))
        one-hot-encoding of DNA sequence where A, C, G, T
        are represented by 0, 1, 2, 3, respectively
    '''
    letter_to_index =  {'A':0, 'a':0,
                        'C':1, 'c':1,
                        'G':2, 'g':2,
                        'T':3, 't':3}
    to_return = np.zeros((len(seq),4), dtype='int8')
    for idx,letter in enumerate(seq):
        if letter in letter_to_index:
            to_return[idx,letter_to_index[letter]] = 1
    return to_return

def encode(seq):
    ''' encode as 0, 1, 2, 3 instead of one-hot-encode for cosine similarity
    '''
    letter_to_index =  {'A':0, 'a':0,
                        'C':1, 'c':1,
                        'G':2, 'g':2,
                        'T':3, 't':3}
    encoding = np.zeros(len(seq))
    for idx, letter in enumerate(seq):
        if letter in letter_to_index:
            encoding[idx] = letter_to_index[letter]
    return encoding


def get_usable_peakids():
  species_peaks = dict()
  peak_list = get_filters("test")
  uniq_species = list()
  for file in tqdm(glob.glob("../../data/200_mammals/*.fa")):
    name = file.split("/")[-1].split("mm10orthlogous_")[1].split("_500bp.fa")[0]
    if name not in uniq_species:
      uniq_species.append(name)
    for seq in SeqIO.parse(file, "fasta"):
      seq_id = seq.id.split(":")[0]
      if seq_id in peak_list:
        if seq_id in species_peaks:
          species_peaks[seq_id].append(name)
        else:
          species_peaks[seq_id] = [name]
  usable_peaks = list()
  for peak in species_peaks:
    if len(species_peaks[peak]) >= 0.25*len(uniq_species):
      usable_peaks.append(peak)
  return usable_peaks


def get_filters(split):
    train_chroms = {"3","4","5","6","7","10","11","12","13","14","15","16","17","18","19","X"}
    val_chroms = {"8", "9"}
    test_chroms = {"1", "2"}
    peak_set = "/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/mouse/GSE146897_WTMouse_ret_ATAC/mm10_ret_noTSS_filtered_labeled.bed"
    #peak_set = "/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/human/GSE137311/hg38_ret_noTSS_filtered_labeled.bed"
    peak_list = list()
    with open(peak_set, "r") as f:
        for line in f.read().splitlines():
            chr_num = line.split("\t")[0].split("chr")[1]
            peak_num = line.split("\t")[3]
            if (split == "val"):
                if chr_num in val_chroms:
                    peak_list.append(peak_num)
            elif (split == "train"):
                if chr_num in train_chroms:
                    peak_list.append(peak_num)
            elif (split == "test"):
                if chr_num in test_chroms:
                    peak_list.append(peak_num)
            else:
                if chr_num in val_chroms or chr_num in train_chroms or chr_num in test_chroms:
                    peak_list.append(peak_num)

    f.close()
    return peak_list

#usable_peaks = get_usable_peakids()

def get_filtered_input(path, mode="val"):
    peak_list = get_filters(mode)
    x = list()
    for seq in SeqIO.parse(path, "fasta"):
      if seq.id.split(":")[0] in usable_peaks:
        x.append(onehot_seq(seq))
        x.append(onehot_seq(seq.reverse_complement()))
    x = np.expand_dims(np.array(x), axis=3)
    print(path, " sequences : " + str(len(x)))
    return x

def get_filtered_seqs(path, mode="val", thresh=0.05):
  peak_list = get_filters(mode)
  out = list()
  total = 0
  for seq in SeqIO.parse(path, "fasta"):
    seq_id = seq.id.split(":")[0]
    count = seq.seq.count("N") + seq.seq.count("n")
    if seq_id in usable_peaks:
      total += 1
      if count <= int(len(seq)*thresh):
        out.append(onehot_seq(seq))
        out.append(onehot_seq(seq.reverse_complement()))
  return np.expand_dims(np.array(out), axis=3), total

def get_filtered_input_encode(path, mode):
    peak_list = get_filters(split=mode)
    x = np.expand_dims(np.array([encode(seq) for seq in SeqIO.parse(path, "fasta") if seq.id[1:].split(":")[0] in peak_list] +
    [encode(seq.reverse_complement()) for seq in SeqIO.parse(path, "fasta") if seq.id.split(":")[0] in peak_list]), axis=2)
    print(path, " sequences : " + str(len(x)))
    return x


def get_input(path, enhancers=None):
    if enhancers:
        x = list()
        ids = list()
        for seq in SeqIO.parse(path, "fasta"):
          peakid = seq.id.split(":")[0]
          if peakid in enhancers:
            x.append(onehot_seq(seq))
            x.append(onehot_seq(seq.reverse_complement()))
            ids.append(peakid)
        x = np.expand_dims(np.array(x), axis=3)
        #print(path, " sequences : " + str(len(x)))
    else:
        x = list()
        ids = list()
        for seq in SeqIO.parse(path, "fasta"):
            peakid = seq.id.split(":")[0]
            x.append(onehot_seq(seq))
            x.append(onehot_seq(seq.reverse_complement()))
            ids.append(peakid)
        x = np.expand_dims(np.array(x), axis=3)
        #print(path, " sequences : " + str(len(x)))
    return x, ids

def score_density_plot(subterr_mammals, model, all_species):
  elife_mammals = np.loadtxt("../../data/elife_mammals.txt", dtype=str)
  fossorial_mammals = np.loadtxt("../../data/fossorial_mammals.txt", dtype=str)
  # store list of overlapping peaks
  overlap_bed = "/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/elife_data/elife-25884-supp10-v1_hg38_overlap.txt"
  overlap_peaks = list()
  with open(overlap_bed, "r") as f:
      overlap_peaks = f.read().splitlines()
  f.close()

  overlap_peaks.append("peak2334") # from overlap with mouse orthologs in human
  overlap_peaks = ["peak13719"]
  #overlap_peaks.append("peak4811")
  #overlap_peaks.append("peak28424")
  # peaks 4811, 28424 come up when overlapping mouse orthologs in human with supp file 11

  with open("mouse4_suppfile10_predictions.txt", "w") as f:
    f.write("Name\t")
    for i in range(len(overlap_peaks)):
      if i < len(overlap_peaks)-1:
        f.write(overlap_peaks[i]+"\t")
      else:
        f.write(overlap_peaks[i]+"\n")

    exclude = ["Chaetophractus_vellerosus", "Dasypus_novemcinctus", "Tolypeutes_matacus", "Ictidomys_tridecemlineatus", "Eptesicus_fuscus", "Myotis_davidii", "Sorex_araneus"]
    maybe_exclude = ["Crocidura_indochinensis", "Elephantulus_edwardii", "Microgale_talazaci", "Mus_pahari", "Tupaia_chinensis", "Tupaia_tana", "Uropsilus_gracilis"]
    all_scores = dict()

    for peak in overlap_peaks:
      subterr_scores = list()
      other_scores = list()
      fossorial_scores = list()
      subterr_species = list()
      names = list()

      for species in tqdm(all_species):
        #if species.split("/")[-1].split("_500bp.fa")[0] in elife_mammals:
        try:
          x, _ = get_input(species, [peak])

          #x = get_filtered_input(species, 'val')
          #scores = model.predict(x).ravel()
          score = sum(model.predict(x).ravel())/2

          name = species.split("/")[-1].split("_500bp.fa")[0]
          if name in all_scores:
            all_scores[name].append(score)
          else:
            all_scores[name] = [score]

          #name = species.split("/")[-1].split("_500bp.fa")[0].split("mm10orthologous_")[1]
          if name in subterr_mammals:
            subterr_species.append(name)
            subterr_scores.append([score])
          elif name not in fossorial_mammals:
            other_scores.append(score)
          elif name in fossorial_mammals:
            fossorial_scores.append(score)
        except:
          name = species.split("/")[-1].split("_500bp.fa")[0]
          if name in all_scores:
            all_scores[name].append("NA")
          else:
            all_scores[name] = ["NA"]
          continue

      if len(subterr_scores) == 0:
        break
      colors = {"Condylura_cristata" : 'green',
                "Chrysochloris_asiatica" : 'red',
                "Heterocephalus_glaber" : 'dodgerblue',
                "Nannospalax_galili" : 'indigo',
                "Scalopus_aquaticus" : 'lightgreen',
                "Echinops_telfairi" : 'lightcoral',
                "Fukomys_damarensis" : 'skyblue',
                "Cricetomys_gambianus" : 'blueviolet',
                "Uropsilus_gracilis" : 'black',
                "Ellobius_lutescens" : 'yellow',
                "Ellobius_talpinus" : 'brown'
                }
      print(all_scores)
      sys.exit()

      other_scores = np.array(other_scores)
      data = {"Predicted Enhancer Activity": other_scores,
              "trait": ["other mammals" for _ in range(len(other_scores))]}
      df = pd.DataFrame(data, columns = ['Predicted Enhancer Activity', 'trait'])
      fig, ax = plt.subplots()
      sns.kdeplot(x="Predicted Enhancer Activity", data=df, ax=ax, color='b')
      data = {"Predicted Enhancer Activity": fossorial_scores,
              "trait": ["fossorial_mammals" for _ in range(len(fossorial_scores))]}
      df = pd.DataFrame(data, columns = ['Predicted Enhancer Activity', 'trait'])
      sns.kdeplot(x="Predicted Enhancer Activity", data=df, ax=ax, color='g')
      for scores, species in zip(subterr_scores, subterr_species):
        ax.scatter(scores, [0.05 for _ in range(len(scores))], label=species, alpha=0.7, color=colors[species])
      # add mean lines
      ymin, ymax = tuple(map(int, ax.get_ylim()))
      mean_other = np.mean(other_scores)
      ax.plot([mean_other for _ in range(ymin, ymax+2)], range(ymin, ymax+2), color='k', linestyle='dashed', label="Mean Predicted Activity of Other Species")
      mean_subterr = np.mean(np.hstack(subterr_scores))
      ax.plot([mean_subterr for _ in range(ymin, ymax+2)], range(ymin, ymax+2), color='b', linestyle='dashed', label="Mean Predicted Activity of Subterranean Species")
      ax.set_xlim(0, 1)
      ax.legend(bbox_to_anchor=(0, -.2), loc='upper left',
          ncol=2, borderaxespad=0)
      plt.savefig("subterrvsall_" + peak + ".png", dpi=300, bbox_inches='tight')
      plt.clf()

    print("writing scores...")
    for name in tqdm(sorted(all_scores)):
      f.write(name + "\t")
      for i in range(len(all_scores[name])):
        if i == (len(all_scores[name]) - 1):
          f.write(str(all_scores[name][i]) + "\n")
        else:
          f.write(str(all_scores[name][i]) + "\t")

def predict_glires_data(model, mode):
  usable_peaks = get_usable_peakids()

  means = dict()
  stds = dict()

  for file in tqdm(glob.glob("../../data/200_mammals/*.fa")):
    name = file.split("/")[-1].split("mm10orthlogous_")[1].split("_500bp.fa")[0]
    #x, _ = get_filtered_seqs(file, mode=mode, thresh=0.05)
    x = get_filtered_input(file, mode=mode)
    scores = model.predict(x).ravel()
    means[name] = np.mean(scores)
    stds[name] = np.std(scores)

  # score mouse
  x, _ = get_input("../../models/mouse4_v3/pos_TEST.fa")
  scores = model.predict(x).ravel()
  means["Mus_musculus"] = np.mean(scores)
  stds["Mus_musculus"] = np.std(scores)

  # load df
  times = dict()
  df = pd.read_excel("../../data/200_mammals_genome_information.xlsx")
  for name in means:
    times[name] = float(df.loc[df['Species'] == name, ['Time Since Split from Mouse (TimeTree)']].to_numpy()[0][0])
  times["Mus_musculus"] = 0.0

  # plot exponential decay figures
  neg_validation, _ = get_input("../../models/mouse4_v3/neg_TEST.fa")
  mean_neg = np.mean(model.predict(neg_validation).ravel())

  def func(x, a, c, d):
    return a*np.exp(-c*x)+d

  popt, pcov = curve_fit(func, [times[name] for name in sorted(times)], [means[name] for name in sorted(times)], p0=(1, 1e-6, 1))
  xx = np.unique(np.array([times[x] for x in times]))
  yy = func(xx, *popt)
  print(yy)

  corr, pval = spearmanr([times[species] for species in sorted(times)], [means[species] for species in sorted(means)])
  print(corr, pval)
  np.save("times", np.array([times[species] for species in sorted(times)]))
  np.save("scores", np.array([stds[species] for species in sorted(stds)]))
  print("saved files")
  # mean decay
  fig = plt.figure()

  plt.plot([times[species] for species in sorted(times)], [means[species] for species in sorted(means)], 'ko', label="Predicted Activity")
  plt.plot([0, int(max(times.values()))], [mean_neg, mean_neg], 'y--', label="Mean Predicted Scores for Negative Validation Set")
  plt.plot(xx, yy, 'r--')
  plt.ylim(top=1)
  plt.ylim(bottom=0)
  plt.xlabel("Divergence from Mouse (MYA)")
  plt.ylabel("Mean of Mouse Enhancer Ortholog Predictions")
  plt.legend()
  plt.savefig("../../plots/glires_"+mode+"_scores.png", dpi=300)
  # std decay
  popt, pcov = curve_fit(func, [times[name] for name in sorted(times)], [stds[name] for name in sorted(stds)], p0=(1, 1e-6, 1))
  xx = np.unique(np.array([times[x] for x in times]))
  yy = func(xx, *popt)
  fig = plt.figure()
  plt.plot([times[species] for species in sorted(times)], [stds[species] for species in sorted(stds)], 'ko', label="Predicted Activity")
  plt.plot(xx, yy, 'r--')
  plt.legend()
  plt.savefig("../../plots/glires_" + mode + "_std_scores.png", dpi=300)


def predict_euarchonta_data(model, mode):


  means = dict()
  stds = dict()

  euarchonta_species = list()
  with open("/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/200_mammals/euarchonta_orthologs.txt", "r") as f:
    data = f.read().splitlines()
    for line in data:
      name = line.split("/")[-1].split("labeled_")[1].split("_summit")[0]
      if name != "Galeopterus_variegatus":
        euarchonta_species.append(name)
  f.close()


  for file in tqdm(glob.glob("/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/human/GSE137311/240mammals/seqs/*.fa")):
    name = file.split("/")[-1].split("_500bp.fa")[0]
    if name in euarchonta_species:
      x = get_filtered_input(file)
      scores = model.predict(x).ravel()
      means[name] = np.mean(scores)
      stds[name] = np.std(scores)

  # score mouse
  x, _ = get_input("../../models/human4/pos_VALIDATION.fa")
  scores = model.predict(x).ravel()
  means["Homo_sapiens"] = np.mean(scores)
  stds["Homo_sapiens"] = np.std(scores)

  # load df
  times = dict()
  df = pd.read_excel("../../data/200_mammals_genome_information.xlsx")
  #for name in means:
  for name in euarchonta_species:
    times[name] = float(df.loc[df['Species'] == name, ['Time Since Split from Human (MYA)']].to_numpy()[0][0])
  times["Homo_sapiens"] = 0.0


  # plot exponential decay figures
  neg_validation, _ = get_input("../../models/human4/neg_VALIDATION.fa")
  mean_neg = np.mean(model.predict(neg_validation).ravel())

  def func(x, a, b):
    return a*np.exp(-b*x)

  popt, pcov = curve_fit(func, [times[name] for name in sorted(times)], [means[name] for name in sorted(times)], p0=(1, 1e-6))
  xx = np.linspace(0, int(max(times.values()))+1, int(max(times.values()))+2)
  yy = func(xx, *popt)
  print(yy)

  # mean decay
  fig = plt.figure()
  plt.plot([times[species] for species in sorted(times)], [means[species] for species in sorted(means)], 'ko', label="Predicted Activity")
  plt.plot([0, int(max(times.values()))], [mean_neg, mean_neg], 'y--', label="Mean Predicted Scores for Negative Validation Set")
  plt.plot(xx, yy, 'r--')
  plt.ylim(top=1)
  plt.ylim(bottom=0)
  plt.xlabel("Divergence from Human (MYA)")
  plt.ylabel("Mean of Human Enhancer Ortholog Predictions")
  plt.legend()
  plt.savefig("../../plots/glires_val_scores.png", dpi=300)
  # std decay
  fig = plt.figure()
  plt.plot([times[species] for species in sorted(times)], [stds[species] for species in sorted(stds)], 'ko')
  plt.savefig("../../plots/glires_val_std_scores.png", dpi=300)

def partition_seqs(peaklist, file):
  overlap, nooverlap = list(), list()
  for seq in SeqIO.parse(file, "fasta"):
    if seq.id.split(":")[0] in peaklist:
      overlap.append(onehot_seq(seq))
      overlap.append(onehot_seq(seq.reverse_complement()))
    else:
      nooverlap.append(onehot_seq(seq))
      nooverlap.append(onehot_seq(seq.reverse_complement()))
  return np.expand_dims(np.array(overlap), axis=3), np.expand_dims(np.array(nooverlap), axis=3)

def make_bar_plot(subterr_mammals, model, mode):

  overlap = dict()
  nooverlap = dict()

  overlap_bed = "/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/elife_data/elife-25884-supp10-v1_hg38_overlap.txt"
  overlap_peaks = list()
  with open(overlap_bed, "r") as f:
      overlap_peaks = f.read().splitlines()
  f.close()
  non_overlap_peaks = list()
  # delete line below for all overlapping peaks
  #overlap_peaks = ['peak13719']
  for seq in SeqIO.parse("../../data/mouse/GSE146897_WTMouse_ret_ATAC/mm10_ret_noTSS_filtered_500bp_labeled.fa", "fasta"):
    peak = seq.id.split(":")[0]
    if peak not in overlap_peaks:
      non_overlap_peaks.append(peak)

  for file in tqdm(glob.glob("../../data/human/GSE137311/240mammals/seqs/*.fa")):
    name = file.split("/")[-1].split("_500bp.fa")[0]
    x_overlap, x_nooverlap = partition_seq(overlap_peaks, file)
    if len(x_overlap) > 0:
      x_overlap = np.array(x_overlap)
      scores_overlap = model.predict(x_overlap).ravel()
      overlap[name] = scores_overlap
    x_nooverlap = np.array(x_nooverlap)
    scores_nooverlap = model.predict(x_nooverlap).ravel()
    nooverlap[name] = scores_nooverlap
    if name in subterr_mammals:
      break

  print(overlap)
  print(nooverlap)
  # score mouse
  #x, ids = get_input("../../data/human/GSE137311/240mammals/seqs/Mus_musculus_500bp.fa")
  #x_overlap = np.array([seq for seq, peakid in zip(x, ids) if peakid in overlap_peaks])
  #x_nooverlap = np.array([seq for seq, peakid in zip(x, ids) if peakid in non_overlap_peaks])
  #scores_overlap = model.predict(x_overlap).ravel()
  #scores_nooverlap = model.predict(x_nooverlap).ravel()
  #overlap["Mus_musculus"] = scores_overlap
  #nooverlap["Mus_musculus"] = scores_nooverlap

  fig = plt.figure()
  #plt.bar(["Accelerated Regions", "Other Regions"], [np.mean(np.hstack(overlap.values())), np.mean(np.hstack(nooverlap.values()))])
  #plt.bar(["OTX2"], [np.mean([overlap[species] for species in overlap if species in subterr_mammals])], label="Subterranean")
  plt.bar(["Accelerated enhancers"], [np.mean([overlap[species] for species in overlap if species not in subterr_mammals])], 'b')
  #plt.errorbar(["OTX2"], [overlap[species] for species in overlap if species not in subterr_mammals])
  plt.bar(["Other enhancers"], [np.mean(np.hstack([nooverlap[species] for species in overlap if species not in subterr_mammals]))], 'r')
  #plt.errorbar(["Other enhancers"], [np.hstack([nooverlap[species] for species in overlap if species not in subterr_mammals])])
  #plt.legend()
  plt.savefig("barplot.png")

def percent_N_evaluation(model, mode):
  threshs = np.linspace(0.0, 0.1, 6)
  percents = list()
  results = {"Fraction of N Threshold":[], "Predicted Activity":[]}
  for thresh in tqdm(threshs):
    seqs, _ = get_filtered_seqs("/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/human/GSE137311/240mammals/seqs/Cavia_aperea_500bp.fa", mode, thresh)
    scores = model.predict(seqs).ravel()
    print(np.mean(scores))
    for score in scores:
      results["Fraction of N Threshold"].append(str(round(thresh, 2)))
      results["Predicted Activity"].append(score)
  df = pd.DataFrame(data=results)
  sns.violinplot(y="Fraction of N Threshold", x="Predicted Activity", data=df, cut=0, palette="Set3")
  plt.savefig("../../plots/threshold_violinplot.png")
  plt.clf()

def count_peaks():
  overlap = np.loadtxt("/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/mouse/GSE146897_WTMouse_ret_ATAC/mouse_orth_overlap_human_peaks.txt", dtype=str)
  count = 0
  
  # for file in tqdm(glob.glob("/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/human/GSE137311/240mammals/seqs/*.fa")):
  #   name = file.split("/")[-1].split("_500bp.fa")[0]
  #   for seq in SeqIO.parse(file, "fasta"):
  #     if (seq.id.split(":")[0]) not in overlap:
  #       count += 1
  for seq in SeqIO.parse("/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/human/GSE137311/240mammals/seqs/Homo_sapiens_500bp.fa", "fasta"):
    if (seq.id.split(":")[0]) not in overlap:
      count += 1
  print(count)


def main(args):

  #MODEL_PATH = "/home/csriniv1/mouse4_v3_mouse4_v3_OCP_NB256_NE10_BR1e-05_MR0.01_BM0.875_MM0.99.h5"
  #MODEL_PATH="/home/csriniv1/mouse4_newarch.h5"
  #MODEL_PATH="/home/csriniv1/mouse4_heinit.h5" # current best model
  #MODEL_PATH="/home/csriniv1/mouse4_heinit2.h5"
  #MODEL_PATH="/home/csriniv1/mouse4_smallarch.h5"
  subterr_mammals = ["Condylura_cristata", "Chrysochloris_asiatica", "Heterocephalus_glaber", "Nannospalax_galili"]
  moles = ["Uropsilus_gracilis", "Ellobius_lutescens", "Ellobius_talpinus", "Fukomys_damarensis", "Scalopus_aquaticus"]
  subterr_neighbors = ["Scalopus_aquaticus", "Echinops_telfairi", "Fukomys_damarensis", "Cricetomys_gambianus"]
  subterr_mammals += subterr_neighbors
  #subterr_mammals += moles
  # "Suricata_suricatta" (meerkat)
  # "Xerus_inauris" (cape ground squirrel)
  # add on fossorial mammals
  #model = load_model(args.model)
  mode = args.mode
  # predict across glires
  count_peaks()

  #predict_glires_data(model, mode)
  #predict_euarchonta_data(model, mode)

  #make_bar_plot(subterr_mammals, model, mode)

  # plot cosine similarity

  #plot_cosine_similarity(mode)
  #cosine_similarity_two_species()

  # plot training/validation score distributions
  #plot_dist(model, species="human", mode="val")
  #plot_dist(model, species="human", mode="train")

  #write_false_negs(model, species="human", mode="val")

  # compare particular species scores

  #human_seqs = "/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/200_mammals/mm10orthologous_Homo_sapiens_500bp.fa"
  #goat_seqs = "/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/200_mammals/mm10orthologous_Capra_hircus_500bp.fa"
  #print(score_seqs(model, human_seqs, "train"))
  #print(score_seqs(model, goat_seqs, "train"))

  # make score density plot

  # score_density_plot(subterr_mammals,
  #                    model,
  #                    #glob.glob("/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/200_mammals/*.fa")
  #                    glob.glob("/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/human/GSE137311/240mammals/seqs/*.fa")
  #                   )
  #percent_N_evaluation(model, mode)



if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='CNN-predicted enhancer scores')
  parser.add_argument("--model", type=str, default="/home/csriniv1/mouse4_heinit2.h5", help="complete model name", required=False)
  parser.add_argument("--mode", type=str, help="Dataset to perform evaluations on.",
        choices=[ 'train', 'val', 'test'], default = 'val', required=False)
  args = parser.parse_args()
  main(args)



