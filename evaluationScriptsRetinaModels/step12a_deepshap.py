from tensorflow.keras.models import load_model
import tensorflow.keras.backend as K
import os, sys, gc, math, argparse
import tensorflow as tf
import pandas as pd
import numpy as np
import shap
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import sys
from sklearn import metrics
from callback_ocp import *
from cnn_utils import *

# gc.collect()
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  # or any {'0', '1', '2'}
# make sure running on GPU #
tf.config.experimental.list_physical_devices('GPU')
# tf.get_logger().setLevel('WARNING')
np.seterr(divide = 'ignore')



def print_2darray(arr):
    s = ''
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            if j == 0 and i == 0: #first entry
                s = str(arr[i,j][0])
            elif j == 0: # new line
                s = s + ';' + str(arr[i,j][0])
            else: # all others
                s = s + ',' + str(arr[i,j][0])
    return s



def print_1darray(arr):
    s = ''
    for i in range(arr.shape[0]):
        s = s + arr[i][0]
    return s



def onehot2fasta(arr,ids,args):
    letters = np.array(['A','C','G','T'])[np.argmax(arr, axis=2)]
    my_seq = [SeqRecord(seq=Seq(print_1darray(letters[i,])), id=ids[i].split('_')[1], description = '',name='') for i in range(len(ids))]
    return my_seq



def predict_sequences(model_name, x, ids):
    # Creating an empty Dataframe with column names only
    # predict labels
    # combineStrands will average the logit of the for and rev DNA strands
    model = load_model(model_name, compile=False)
    y_pred_score = model.predict(x).flatten()
    y_pred_score = np.nan_to_num(y_pred_score)
    df = pd.DataFrame({'id': ids, 'y_pred_logit': np.log10(y_pred_score) - np.log10(1-y_pred_score)})
    df = df.groupby(['id']).mean()
    df['y_pred_score'] = 1 / (1 + np.exp(-df['y_pred_logit']))
    df['y_pred_class'] = df['y_pred_score'] > 0.5
    return df



def evaluate_SHAP_scores(args, fg, fg_ids, bg):
    ## create deep SHAP explainer for bg vs. fg
    minRepBg = np.min([bg.shape[0], args.numBg])
    minRepFg = np.min([fg.shape[0], args.numFg])
    np.random.seed(seed=args.seed)
    bg = bg[np.random.choice(bg.shape[0], minRepBg, replace = False)]
    indFg = np.random.choice(fg.shape[0], minRepFg, replace = False)
    fg = fg[indFg]; ids = fg_ids[indFg]
    chunks = 10
    fg_list = np.array_split(fg, fg.shape[0] // chunks)
    ## initialize deep SHAP explainer for bg vs. fg & compute fg in 100 seq chunks
    model = load_model(args.model_name, compile=False)
    explainer = shap.DeepExplainer(model, bg)
    rawShapList = []
    for ind in range(0,len(fg_list)):
        print(f'Computing SHAP values: {(ind + 1) * chunks} of {fg.shape[0]}.')
        raw_shap_explanations = explainer.shap_values(fg_list[ind], check_additivity = False)
        rawShapList.append(raw_shap_explanations)
    # combine explanations and normalize
    hyp = np.squeeze(np.concatenate( rawShapList, axis = 1 ), axis=0)
    imp = hyp*fg # mask hypothetical w/ actual sequences
    return (ids, hyp, imp, fg)



def main(args):
    """Main function
    Args:
    args (argparse):
    """
    # file names
    # call main functions
    if not os.path.exists(args.model_name):
        print('No model found with specified training parameters. Please train model first.')
        return
    if not os.path.exists(f'{args.out_dir}/shap'):
        os.makedirs(f'{args.out_dir}/shap')
    if args.mode == 'evaluate':
        print('In evaluation mode to compute SHAP scores.')
        print(f'Reading in positive sequences from {args.eval_fasta_pos}.')
        (x_pos, ids_pos) = encode_sequence3(args.eval_fasta_pos, size = args.seq_length, shuffleOff = True)
        x_pos = x_pos[range(len(ids_pos)//2),:]; ids_pos = ids_pos[range(len(ids_pos)//2)]
        df_pos = predict_sequences(args.model_name, x_pos, ids_pos).loc[ids_pos,:]

        indKeep_pos = np.where(df_pos['y_pred_class'] == 1)[0]
        

        #x_pos = x_pos[indKeep_pos]; ids_pos = ids_pos[indKeep_pos
        # uncomment for only positive seqs 
        df_pos = df_pos.loc[ids_pos,:]
        #
        print(f'Reading in negative sequences from {args.eval_fasta_neg}.')
        (x_neg, ids_neg) = encode_sequence3(args.eval_fasta_neg, size = args.seq_length, shuffleOff = True)
        x_neg = x_neg[range(len(ids_neg)//2),:]; ids_neg = ids_neg[range(len(ids_neg)//2)]
        df_neg = predict_sequences(args.model_name, x_neg, ids_neg).loc[ids_neg,:]
        indKeep_neg = np.where(df_neg['y_pred_class'] == 0)[0]
        x_neg = x_neg[indKeep_neg]; ids_neg = ids_neg[indKeep_neg]
        df_neg = df_neg.loc[ids_neg,:]
        #
        # get hypothetical and importance scores
        print(f'Computing SHAP scores with true positives and true negatives.')
        (ids, hyp, imp, fg) = evaluate_SHAP_scores(args, fg = x_pos, fg_ids = ids_pos, bg = x_neg)
        df = df_pos.loc[ids,:]
        #
        # save hypothetical importance scores
        imp_out = pd.DataFrame({'id': df.index.str.split(pat = '_').str[1],
            'y_pred_logit': df['y_pred_logit'],
            'SHAP_score': [print_2darray(imp[i,:]) for i in range(imp.shape[0])]})
        imp_fn = f'{args.out_dir}/shap/{args.predict_out}.pos{args.numFg}.imp_SHAP_scores.txt'
        print(f'Writing importance scores to {imp_fn}.')
        imp_out.to_csv(imp_fn, index = False, sep = '\t', header = False)
        #
        hyp_out = pd.DataFrame({'id': df.index.str.split(pat = '_').str[1],
            'y_pred_logit': df['y_pred_logit'],
            'SHAP_score': [print_2darray(hyp[i,:])  for i in range(hyp.shape[0])]})
        hyp_fn = f'{args.out_dir}/shap/{args.predict_out}.pos{args.numFg}.hyp_SHAP_scores.txt'
        print(f'Writing hypothetical importance scores to {hyp_fn}.')
        hyp_out.to_csv(hyp_fn, index = False, sep = '\t', header = False)
        #
        fasta_out = onehot2fasta(fg,ids,args)
        fasta_fn = f'{args.out_dir}/shap/{args.predict_out}.pos{args.numFg}.fasta'
        print(f'Writing scored DNA sequences {fasta_fn}.')
        SeqIO.write(fasta_out, fasta_fn, "fasta")
        #
    elif args.mode == 'predict':
        print('In prediction mode to compute SHAP scores.')
        print(f'Reading in sequences from {args.predict_fasta}.')
        (x, ids) = encode_sequence3(args.predict_fasta, size = args.seq_length, shuffleOff = True)
        df = predict_sequences(args.model_name, x, ids).loc[ids,:]
        ids = np.array([ ids[i] if i < len(ids)//2 else ids[i] + '.rev' for i in range(len(ids)) ])
        df.index = ids; args.numFg = len(ids)
        #
        print(f'Reading in negative sequences from {args.eval_fasta_neg}.')
        (x_neg, ids_neg) = encode_sequence3(args.eval_fasta_neg, size = args.seq_length, shuffleOff = True)
        x_neg = x_neg[range(len(ids_neg)//2),:]; ids_neg = ids_neg[range(len(ids_neg)//2)]
        df_neg = predict_sequences(args.model_name, x_neg, ids_neg).loc[ids_neg,:]
        indKeep_neg = np.where(df_neg['y_pred_class'] == 0)[0]
        x_neg = x_neg[indKeep_neg]; ids_neg = ids_neg[indKeep_neg]
        df_neg = df_neg.loc[ids_neg,:]
        #
        # get hypothetical and importance scores
        print(f'Computing SHAP scores of sequences against true negatives.')
        (ids, hyp, imp, fg) = evaluate_SHAP_scores(args, fg = x, fg_ids = ids, bg = x_neg)
        df = df.loc[ids,:]
        #
        # save hypothetical importance scores
        imp_out = pd.DataFrame({'id': df.index.str.split(pat = '_').str[1],
            'y_pred_logit': df['y_pred_logit'],
            'SHAP_score': [print_2darray(imp[i,:]) for i in range(imp.shape[0])]})
        imp_fn = f'{args.out_dir}/shap/{args.predict_out}.predict.imp_SHAP_scores.txt'
        print(f'Writing importance scores to {imp_fn}.')
        imp_out.to_csv(imp_fn, index = False, sep = '\t', header = False)
        #
        hyp_out = pd.DataFrame({'id': df.index.str.split(pat = '_').str[1],
            'y_pred_logit': df['y_pred_logit'],
            'SHAP_score': [print_2darray(hyp[i,:])  for i in range(hyp.shape[0])]})
        hyp_fn = f'{args.out_dir}/shap/{args.predict_out}.predict.hyp_SHAP_scores.txt'
        print(f'Writing hypothetical importance scores to {hyp_fn}.')
        hyp_out.to_csv(hyp_fn, index = False, sep = '\t', header = False)
        #
        fasta_out = onehot2fasta(fg,ids,args)
        fasta_fn = f'{args.out_dir}/shap/{args.predict_out}.pos{args.numFg}.fasta'
        print(f'Writing scored DNA sequences {fasta_fn}.')
        SeqIO.write(fasta_out, fasta_fn, "fasta")
    return



if __name__ == '__main__':
    #### set cnn parameters:
    parser = argparse.ArgumentParser(description='Parse CNN parameters.')
    parser.add_argument("--mode", type=str, help="Mode to perform. Train needs all fasta. Evaluate needs validation fasta. Predict only fasta passed predict_fasta.",
        choices=[ 'evaluate', 'predict'], default = 'evaluate', required=False)
    #
    parser.add_argument("--model_name", type=str, help="complete model name", required=True)
    parser.add_argument("--predict_out", type=str, help="prediction file prefix model name", required=True)
    parser.add_argument("--predict_fasta", type=str, help="fasta sequence file for predictions.", required=False)
    parser.add_argument("--eval_fasta_pos", type=str, help="validation fasta sequence file of positives.", required=False)
    parser.add_argument("--eval_fasta_neg", type=str, help="validation fasta sequence file of negatives.", required=True)
    parser.add_argument("--numFg", type=int, default = 2000, help="number of positives to evaluate with SHAP scores. Note that all predict fasta are scored.", required=False)
    parser.add_argument("--numBg", type=int, default = 500, help="number of negatives in background to compute SHAP scores", required=False)
    parser.add_argument("--seq_length", type=str, default = 500, help="number of negatives in background to compute SHAP scores", required=False)
    parser.add_argument("--seed", type=int, default = 1, help="number of negatives", required=False)
    parser.add_argument("--force", help="Whether to overwrite previously trained model.", action='store_true')
    parser.add_argument("--verbose", type=int, default = 2, help="Verbosity in keras.")
    parser.add_argument("--out_dir", type=str, default = '.', help="path to ouputput directory, default is pwd")

    # args = parser.parse_args(['--model_name=/projects/pfenninggroup/machineLearningForComputationalBiology/addiction_gwas_enrichment/celltype_specific_ml_models/models/Mo2015_EXCpos_Ctx_fold1/Mo2015_EXCpos_Ctx_fold1_OCP_NB1000_NE23_BR0.01_MR0.1_BM0.85_MM0.99_DO0.25.h5',
    #     '--predict_out=Mo2015_EXCpos_Ctx_fold1',
    #     '--eval_fasta_pos=FASTA_CV/Mo2015_EXCpos_Ctx_fold1_validPos.fa',
    #     '--eval_fasta_neg=FASTA_CV/Mo2015_EXCpos_Ctx_fold1_validNeg10x.fa',
    #     '--mode=evaluate'])

    # args = parser.parse_args(['--model_name=/projects/pfenninggroup/machineLearningForComputationalBiology/addiction_gwas_enrichment/celltype_specific_ml_models/models/Mo2015_EXCpos_Ctx_fold1/Mo2015_EXCpos_Ctx_fold1_OCP_NB1000_NE23_BR0.01_MR0.1_BM0.85_MM0.99_DO0.25.h5',
    #     '--predict_out=top_addiction_snps_effect_allele.Mo2015_EXCpos_Ctx_fold1',
    #     '--predict_fasta=FASTA/top_addiction_snps_allele_effect_501.fa',
    #     '--eval_fasta_neg=FASTA_CV/Mo2015_EXCpos_Ctx_fold1_validNeg10x.fa',
    #     '--mode=predict'])

    # args.numBg = 100

    args = parser.parse_args()
    main(args)



'''
command to run 06/15/21
python step12a_deepshap.py --model_name /projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/mouse4_v3/tmp.h5 --predict_out test2 --eval_fasta_pos /projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/mouse4_v3/pos_VALIDATION.fa --eval_fasta_neg /projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/mouse4_v3/neg_VALIDATION.fa
'''
