#runModisco.py
# usage: runModisco.py <fastaFile> <impScoreFile> <hypScoreFile> runName

import modisco
from matplotlib import pyplot as plt
plt.style.use('default')
import numpy as np
import vizsequence
import sys, os
import h5py

def one_hot_encode_along_channel_axis(sequence):
    to_return = np.zeros((len(sequence),4), dtype=np.int8)
    seq_to_one_hot_fill_in_array(zeros_array=to_return,
                                 sequence=sequence, one_hot_axis=1)
    return to_return



def seq_to_one_hot_fill_in_array(zeros_array, sequence, one_hot_axis):
    assert one_hot_axis==0 or one_hot_axis==1
    if (one_hot_axis==0):
        assert zeros_array.shape[1] == len(sequence)
    elif (one_hot_axis==1):
        assert zeros_array.shape[0] == len(sequence)
    #will mutate zeros_array
    for (i,char) in enumerate(sequence):
        if (char=="A" or char=="a"):
            char_idx = 0
        elif (char=="C" or char=="c"):
            char_idx = 1
        elif (char=="G" or char=="g"):
            char_idx = 2
        elif (char=="T" or char=="t"):
            char_idx = 3
        elif (char=="N" or char=="n"):
            continue #leave that pos as all 0's
        else:
            raise RuntimeError("Unsupported character: "+str(char))
        if (one_hot_axis==0):
            zeros_array[char_idx,i] = 1
        elif (one_hot_axis==1):
            zeros_array[i,char_idx] = 1



def normalize_scores(impscores, hyp_impscores, onehotData):
  normed_hyp_impscores = []
  normed_impscores = []
  for i in range(len(impscores)):
      imp_score_each_pos = np.sum(impscores[i],axis=-1)
      imp_score_sign_each_pos = np.sign(imp_score_each_pos)
      hyp_scores_same_sign_mask = (np.sign(hyp_impscores[i])
                                   *imp_score_sign_each_pos[:,None] > 0)
      hyp_scores_same_sign_imp_scores_sum = np.sum(
          hyp_impscores[i]*hyp_scores_same_sign_mask,axis=-1)
      norm_ratio = imp_score_each_pos/hyp_scores_same_sign_imp_scores_sum
      norm_hyp = hyp_impscores[i]*norm_ratio[:,None]
      normed_hyp_impscores.append(norm_hyp)
      normed_impscores.append(norm_hyp*onehotData[i])
  return normed_impscores, normed_hyp_impscores


def main():
  fastaFile = sys.argv[1]
  impFile = sys.argv[2]
  hypFile = sys.argv[3]
  runName = sys.argv[4]

  os.system("mkdir -p pwm_full")
  os.system("mkdir -p pwm_crop11")
  os.system("mkdir -p results")
  fasta_seqs = [x.rstrip() for (i,x) in enumerate(open(fastaFile))
              if i%2==1]
  onehot_data = np.array([one_hot_encode_along_channel_axis(x)
                           for x in fasta_seqs])
  impscores = [
       np.array( [[float(z) for z in y.split(",")] for y in x.rstrip().split("\t")[2].split(";")])
       for x in open(impFile)
  ]

  hyp_impscores = [
       np.array( [[float(z) for z in y.split(",")] for y in x.rstrip().split("\t")[2].split(";")])
       for x in open(hypFile)
  ]

  print(len(impscores))
  print(len(onehot_data))
  print(len(onehot_data[1]))

  normed_impscores, normed_hyp_impscores = normalize_scores(
  impscores=impscores, hyp_impscores=hyp_impscores, onehotData=onehot_data)

  tfmodisco_results = modisco.tfmodisco_workflow.workflow.TfModiscoWorkflow(
                           sliding_window_size=11,
                           flank_size=3,
                           min_seqlets_per_task=3000,
                           min_passing_windows_frac=0.005, #default
                           seqlets_to_patterns_factory=modisco.tfmodisco_workflow.seqlets_to_patterns.TfModiscoSeqletsToPatternsFactory(
                           trim_to_window_size=11,
                           initial_flank_to_add=3,
                           final_flank_to_add=3,
                           kmer_len=7, num_gaps=1,
                           num_mismatches=1,
                           ),
                      )(
                   task_names=["task0"],
                   contrib_scores={'task0': normed_impscores},
                   hypothetical_contribs={'task0': normed_hyp_impscores},
                   one_hot=onehot_data)

  o_pwm = "pwm_full/%s_pwm_full.txt" % runName
  o_h5 = "results/%s_tfmodiscoResults.hdf5" % runName
  grp = h5py.File(o_h5, 'w')
  tfmodisco_results.save_hdf5(grp)
  grp.close()


  with open(o_pwm, 'w') as f:
      f.write("MEME version 4" + "\n")
      f.write("\n")
      f.write("ALPHABET= ACGT" + "\n")
      f.write("\n")
      f.write("strands: + -" + "\n")
      f.write("\n")
      f.write("Background letter frequencies (from unknown source)" + "\n")
      f.write("A 0.250 C 0.250 G 0.250 T 0.250" + "\n")
      f.write("\n")

      hdf5_results = h5py.File(o_h5,"r")
      metacluster_names = [
          x.decode("utf-8") for x in
          list(hdf5_results["metaclustering_results"]
               ["all_metacluster_names"][:])]
      background = np.mean(onehot_data, axis=(0,1))
      for metacluster_name in metacluster_names:
        print(metacluster_name)
        metacluster_grp = (hdf5_results["metacluster_idx_to_submetacluster_results"]
                                       [metacluster_name])
        print("activity pattern:",metacluster_grp["activity_pattern"][:])
        all_pattern_names = [x.decode("utf-8") for x in
                             list(metacluster_grp["seqlets_to_patterns_result"]
                                                 ["patterns"]["all_pattern_names"][:])]
        print(all_pattern_names)
        print(hdf5_results.metacluster_idx_to_submetacluster_results[0].seqlets_to_patterns_result.patterns)
        sys.exit()
        for i, pattern in enumerate(all_pattern_names):
          #if pattern in enumerate(tfmodisco_results.metacluster_idx_to_submetacluster_results[0].seqlets_to_patterns_result.patterns):
          f.write("MOTIF " + runName + "_"+ str(i+1) + ".fwd" "\n")
          try:
            f.write("letter-probability matrix: alength= 4 w= 21 nsites= " + str(len(pattern.seqlets)) + " E= 0\n")
            for b in range(0,21):
              f.write(str(pattern["sequence"].fwd[b][0]) + "\t" + str(pattern["sequence"].fwd[b][1]) + "\t" + str(pattern["sequence"].fwd[b][2]) + "\t" + str(pattern["sequence"].fwd[b][3]) + "\n")
          except:
            f.write("letter-probability matrix: alength= 4 w= 21 nsites= " + str(0) + " E= 0\n")
        f.write("\n")


  o_pwm = "pwm_full/%s_pwm_full.txt" % runName
  with open(o_pwm, 'w') as f:
      f.write("MEME version 4" + "\n")
      f.write("\n")
      f.write("ALPHABET= ACGT" + "\n")
      f.write("\n")
      f.write("strands: + -" + "\n")
      f.write("\n")
      f.write("Background letter frequencies (from unknown source)" + "\n")
      f.write("A 0.250 C 0.250 G 0.250 T 0.250" + "\n")
      f.write("\n")
      for i,pattern in enumerate(tfmodisco_results.metacluster_idx_to_submetacluster_results[0].seqlets_to_patterns_result.patterns):
          f.write("MOTIF " + runName + "_"+ str(i+1) + ".fwd" "\n")
          f.write("letter-probability matrix: alength= 4 w= 21 nsites= " + str(len(pattern.seqlets)) + " E= 0")
          f.write("\n")
          for b in range(0,21):
              f.write(str(pattern["sequence"].fwd[b][0]) + "\t" + str(pattern["sequence"].fwd[b][1]) + "\t" + str(pattern["sequence"].fwd[b][2]) + "\t" + str(pattern["sequence"].fwd[b][3]) + "\n")
          f.write("\n")



  o_pwmc = "pwm_crop11/%s_pwm_crop11.txt" % runName
  with open(o_pwmc, 'w') as f:
      f.write("MEME version 4" + "\n")
      f.write("\n")
      f.write("ALPHABET= ACGT" + "\n")
      f.write("\n")
      f.write("strands: + -" + "\n")
      f.write("\n")
      f.write("Background letter frequencies (from unknown source)" + "\n")
      f.write("A 0.250 C 0.250 G 0.250 T 0.250" + "\n")
      f.write("\n")
      for i,pattern in enumerate(tfmodisco_results.metacluster_idx_to_submetacluster_results[0].seqlets_to_patterns_result.patterns):
          f.write("MOTIF " + runName + "_" + str(i+1) + ".fwd" "\n")
          f.write("letter-probability matrix: alength= 4 w= 11 nsites= " + str(len(pattern.seqlets)) + " E= 0")
          f.write("\n")
          for b in range(5,16):
              f.write(str(pattern["sequence"].fwd[b][0]) + "\t" + str(pattern["sequence"].fwd[b][1]) + "\t" + str(pattern["sequence"].fwd[b][2]) + "\t" + str(pattern["sequence"].fwd[b][3]) + "\n")
          f.write("\n")

  o_h5 = "results/%s_tfmodiscoResults.hdf5" % runName
  grp = h5py.File(o_h5, 'w')
  tfmodisco_results.save_hdf5(grp)
  grp.close()

  print("done")


  '''
  example command:
  python step12b_modisco.py shap/test2.pos2000.fasta shap/test2.pos2000.imp_SHAP_scores.txt shap/test2.pos2000.hyp_SHAP_scores.txt test2
  '''

if __name__ == "__main__":
    main()


