'''
Environment (Python 3.6.12)
- Keras (2.3.0-tf)
- TensorFlow (2.2.0)
- Numpy (1.19.2)
- Bio (1.78)
'''

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

from sklearn.metrics import auc, roc_curve, precision_recall_curve
import matplotlib.pyplot as plt

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
        if letter not in ['N','n']:
            to_return[idx,letter_to_index[letter]] = 1
    return to_return

def write_false_negs(idxs, infile, outfile1, outfile2):
    mouse_data="/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/mouse/GSE146897_WTMouse_ret_ATAC/"
    human_data="/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/human/GSE137311/"
    eval1 = mouse_data+"mouse_specific_mm10.fa"
    eval2 = mouse_data+"mouse_specific_hg38.fa"
    eval3 = human_data+"human_specific_hg38.fa"
    eval4 = human_data+"human_specific_mm10.fa"
    eval5 = mouse_data+"mm10_ret_nooverlap_liver.fa"
    eval6 = mouse_data+"mm10_ret_overlap_liver.fa"
    eval7 = mouse_data+"mm10_ret_noTSS_filtered_500bp_neg_overlap_liver.fa"
    with open(eval2, "r") as f:
        data = f.read().splitlines()
    f.close()
    counter = 0
    with open(outfile1, "w") as f, open(outfile2, "w") as g:
        for i in tqdm(range(0, len(data), 2)):
            if counter in idxs:
                f.write(data[i] +"\n")
                f.write(data[i+1]+"\n")
            else:
                g.write(data[i] +"\n")
                g.write(data[i+1]+"\n")
            counter += 1
    f.close()
    g.close()

def get_input(path):
    x = np.expand_dims(np.array([onehot_seq(seq) for seq in SeqIO.parse(path, "fasta") ] +
    [onehot_seq(seq.reverse_complement()) for seq in SeqIO.parse(path, "fasta") ]), axis=3)
    y = np.ones(len(x))
    return x, y

def get_filtered_input(path):
    val_chroms = {"chr8", "chr9"}
    x = np.expand_dims(np.array([onehot_seq(seq) for seq in SeqIO.parse(path, "fasta") if seq.id.split(":")[0] in val_chroms] +
    [onehot_seq(seq.reverse_complement()) for seq in SeqIO.parse(path, "fasta") if seq.id.split(":")[0] in val_chroms]), axis=3)
    y = np.ones(len(x))
    print(path, " sequences : " + str(len(y)))
    return x, y

def calc_acc(metrics):
    return (metrics[1] + metrics[3])/sum(metrics[1:5])

def load_data(filter, model):
    mouse_data="/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/mouse/GSE146897_WTMouse_ret_ATAC/"
    human_data="/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/human/GSE137311/"
    eval1 = mouse_data+"mouse_specific_mm10.fa"
    eval2 = mouse_data+"mouse_specific_hg38.fa"
    eval3 = human_data+"human_specific_hg38.fa"
    eval4 = human_data+"human_specific_mm10.fa"
    eval5 = mouse_data+"mm10_ret_nooverlap_liver.fa"
    eval6 = mouse_data+"mm10_ret_overlap_liver.fa"
    eval7 = mouse_data+"mm10_ret_noTSS_filtered_500bp_neg_overlap_liver.fa"
    if filter:
        x_eval1, y_eval1 = get_filtered_input(eval1)
        x_eval2, _ = get_filtered_input(eval2)
        y_eval2 = np.zeros(len(x_eval2))
        x_eval3, y_eval3 = get_filtered_input(eval3)
        x_eval4, _ = get_filtered_input(eval4)
        y_eval4 = np.zeros(len(x_eval4))
        x_eval5, y_eval5 = get_filtered_input(eval5)
        x_eval6, y_eval6 = get_filtered_input(eval6)
        x_eval7, _ = get_filtered_input(eval7)
        y_eval7 = np.zeros(len(x_eval7))
    else:
        x_eval1, y_eval1 = get_input(eval1)
        x_eval2, _ = get_input(eval2)
        y_eval2 = np.zeros(len(x_eval2))
        x_eval3, y_eval3 = get_input(eval3)
        x_eval4, _ = get_input(eval4)
        y_eval4 = np.zeros(len(x_eval4))
        x_eval5, y_eval5 = get_input(eval5)
        x_eval6, y_eval6 = get_input(eval6)
        x_eval7, _ = get_input(eval7)
        y_eval7 = np.zeros(len(x_eval7))
    metrics1 = model.evaluate(x_eval1, y_eval1)
    metrics2 = model.evaluate(x_eval2, y_eval2)
    metrics3 = model.evaluate(x_eval3, y_eval3)
    metrics4 = model.evaluate(x_eval4, y_eval4)
    metrics5 = model.evaluate(x_eval5, y_eval5)
    metrics6 = model.evaluate(x_eval6, y_eval6)
    metrics7 = model.evaluate(x_eval7, y_eval7)
    scores2 = model.predict(x_eval2).ravel()
    scores3 = model.predict(x_eval3).ravel()
    false_negs2 = write_false_negs(np.argwhere(scores2 >= 0.5),
                                    eval2, eval2[:-3]+"_fn.fa", eval2[:-3]+"_tp.fa")
    print(eval2[:-3]+"_fn.fa")
    false_negs3 = write_false_negs(np.argwhere(scores3 < 0.5),
                                    eval3, eval3[:-3]+"_fn.fa", eval3[:-3]+"_tp.fa")
    print(eval3[:-3]+"_fn.fa")
    print(calc_acc(metrics1))
    print(calc_acc(metrics2))
    print(calc_acc(metrics3))
    print(calc_acc(metrics4))
    print(calc_acc(metrics5))
    print(calc_acc(metrics6))
    print(calc_acc(metrics7))
    return [[np.concatenate((x_eval1,x_eval4)), np.concatenate((y_eval1,y_eval4))],
            [np.concatenate((x_eval3,x_eval2)), np.concatenate((y_eval3,y_eval2))],
            #[np.concatenate((x_eval1,x_eval2)), np.concatenate((y_eval1,y_eval2))],
            #[np.concatenate((x_eval3,x_eval4)), np.concatenate((y_eval3,y_eval4))],
            [np.concatenate((x_eval5,x_eval7)), np.concatenate((y_eval5,y_eval7))],
            [np.concatenate((x_eval6,x_eval7)), np.concatenate((y_eval6,y_eval7))]]

def make_roc(model, x, y, name):
    fig = plt.figure()
    plt.plot([0, 1], [0, 1], 'k--')
    y_pred = model.predict(x).ravel()
    fpr_keras, tpr_keras, thresholds_keras = roc_curve(y, y_pred)
    auc_keras = auc(fpr_keras, tpr_keras)
    plt.plot(fpr_keras, tpr_keras, label='(eval area = {:.3f})'.format(auc_keras))
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.title(name + ' ROC curve')
    plt.legend(loc='best')
    plt.savefig("evaluation_curves/" + name + "_roc.png")
    return "evaluation_curves/" + name + "_roc.png"

def make_prc(model, x, y, name):
    no_skill = len(y[y==0]) / len(y)
    fig = plt.figure()
    plt.plot([0, 1], [no_skill, no_skill], linestyle='--', label='No Skill')
    y_pred = model.predict(x).ravel()
    precision, recall, thresholds = precision_recall_curve(y, y_pred)
    auc_keras = auc(recall, precision)
    plt.plot(recall, precision, marker='.', label='(eval area = {:.3f})'.format(auc_keras))
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title(name + ' PRC curve')
    plt.legend(loc='best')
    plt.savefig("evaluation_curves/" + name + "_prc.png")
    return "evaluation_curves/" + name + "_prc.png"

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

def make_plots(model, evals):
    names = ["1_4", "2_3", "1_2", "3_4", "5_7", "6_7"]
    for i in range(len(evals)):
        path = make_roc(model,
                evals[i][0],
                evals[i][1],
                names[i])
        print(path)
        path = make_prc(model,
                evals[i][0],
                evals[i][1],
                names[i])
        print(path)

def bar_plots(model, evals):
    fig = plt.subplot(111)
    names = ["Mouse-specific enhancers", "Human-specific enhancers", "Retina-specific enhancers", "Tissue-shared enhancers"]
    for i in range(len(evals)):
        x, y = evals[i][0], evals[i][1]
        y_pred = model.predict(x).ravel()
        precision, recall, thresholds = precision_recall_curve(y, y_pred)
        auprc = auc(recall, precision)
        fpr_keras, tpr_keras, thresholds_keras = roc_curve(y, y_pred)
        auroc = auc(fpr_keras, tpr_keras)
        plt.bar([i*0.8, 5+(i*0.8)], [auprc, auroc], label=names[i])
    #tick_labels = ['0','AUPRC','5','AUROC']
    #plt.xticks([0, 5], [''.join(x) for x in zip(tick_labels[0::2], tick_labels[1::2])])
    plt.tick_params(bottom=False, labelbottom=False)
    #plt.xlabel("AUPRC" + " "*5 + "AUROC", loc='left')
    plt.ylabel("AUC")
    plt.ylim(0, 1)
    plt.legend(bbox_to_anchor=(0, -.1), loc='upper left',
          ncol=2, borderaxespad=0)
    plt.savefig("test.png", dpi=300, bbox_inches='tight')

def main(options):
    print("loading model...")
    model_path="/home/csriniv1/mouse4_v3_mouse4_v3_OCP_NB256_NE10_BR1e-05_MR0.01_BM0.875_MM0.99.h5"
    model_path="/home/csriniv1/mouse4_2layers.h5"
    model_path="/home/csriniv1/mouse4_heinit.h5"
    model_path="/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/mouse4_v3/tmp.h5"
    if options.model:
        model_path = options.model
    model = load_model(model_path, custom_objects={"macro_f1": macro_f1})

    print("loading sequences...")
    evals = load_data(False, model)
    if options.plot:
        #make_plots(model, evals)
        bar_plots(model, evals)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--model", type=str, help="Path to trained TensorFlow model")
    parser.add_argument("-p", "--plot", type=int, help="Path to plot")

    options, args = parser.parse_known_args()
    '''
    if (options.write):
        parser.print_help(sys.stderr)
        sys.exit(1)
    '''
    main(options)
