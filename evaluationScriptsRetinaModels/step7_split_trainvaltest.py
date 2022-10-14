import argparse
import os
import sys
import numpy as np

# split human, mouse retina ATAC-seq 500 bp summit-centered sequences into train, validation, test sets
# split human data based on HALPER chromosome

# note: negative FASTA sequences obtained from biasaway have subsequences of the 500 bp on newlines

halper_file = "/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/human/GSE137311/hg38_ret_noTSS_filtered_halper_mm10.bed"

human_pos_bed = "/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/human/GSE137311/hg38_ret_noTSS_filtered_labeled_500bp.bed"
human_pos = "/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/human/GSE137311/hg38_ret_noTSS_filtered_labeled_500bp.fa"
human_neg = "/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/human/GSE137311/hg38_ret_noTSS_filtered_500bp_neg_nopos.fa"

human_data="/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/human/GSE137311/"
mouse_data="/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/mouse/GSE146897_WTMouse_ret_ATAC/"
eval4 = human_data+"human_specific_mm10.fa"
eval2 = mouse_data+"mouse_specific_hg38.fa"

mouse_pos ="/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/mouse/GSE146897_WTMouse_ret_ATAC/mm10_ret_noTSS_filtered_500bp.fa"
mouse_neg ="/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/mouse/GSE146897_WTMouse_ret_ATAC/mouse_neg.fa"
mouse_neg_cortex = "/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/mouse/GSE146897_WTMouse_ret_ATAC/mm10_ret_noTSS_filtered_500bp_neg_overlap_cortex.fa"
mouse_human_spec = "/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/human/GSE137311/human_specific_mm10.fa"

train_chroms = {"3","4","5","6","7","10","11","12","13","14","15","16","17","18","19","X"}
val_chroms = {"8", "9"}
test_chroms = {"1", "2"}

# find human dataset splits based on HALPER
def make_human_splits(halper_file, model_path, ratio):
  print("splitting human data...")
  # assign unique orthologous peak identifiers to dataset
  train_set, val_set, test_set = set(), set(), set()
  with open(halper_file, "r") as f:
    contents = f.read().splitlines()
  f.close()
  for line in contents:
    data = line.split("\t")
    if data[0][-1] in train_chroms:
      train_set.add(data[4])
    elif data[0][-1] in val_chroms:
      val_set.add(data[4])
    elif data[0][-1] in test_chroms:
      test_set.add(data[4])

  # match unique identifiers to BED coordinates
  with open(human_pos_bed, "r") as f:
    contents = f.read().splitlines()
  f.close()
  train_peak, val_peak, test_peak = set(), set(), set()
  for line in contents:
    data = line.split("\t")
    if data[3] in train_set:
      # fasta format for matching
      train_peak.add(">"+data[0]+":"+data[1]+"-"+data[2])
    elif data[3] in val_set:
      val_peak.add(">"+data[0]+":"+data[1]+"-"+data[2])
    elif data[3] in test_set:
      test_peak.add(">"+data[0]+":"+data[1]+"-"+data[2])

  # match BED coordinates to FASTA
  with open(human_pos, "r") as f:
    contents = f.read().splitlines()
  f.close()
  with open(human_neg, "r") as f:
    ncontents = f.read().splitlines()
  f.close()
  with open(eval2, "r") as f:
    eval4contents = f.read().splitlines()
  f.close()
  with open(mouse_neg_cortex, "r") as f:
    ncortex = f.read().splitlines()
  f.close()


  with open(model_path+"/pos_TRAINING.fa", "w") as ptr, open(model_path+"/neg_TRAINING.fa", "w") as ntr,   open(model_path+"/pos_VALIDATION.fa", "w") as pva, open(model_path+"/neg_VALIDATION.fa", "w") as nva, open(model_path+"/pos_TEST.fa", "w") as pte, open(model_path+"/neg_TEST.fa", "w") as nte:
    # positive set
    tr_count, va_count, te_count = 0, 0, 0
    for i in range(0, len(contents), 2):
      if "N" not in contents[i+1].upper():
        if contents[i] in train_peak:
          ptr.write(contents[i]+"\n")
          ptr.write(contents[i+1]+"\n")
          tr_count += 1
        elif contents[i] in val_peak:
          pva.write(contents[i]+"\n")
          pva.write(contents[i+1]+"\n")
          va_count += 1
        elif contents[i] in test_peak:
          pte.write(contents[i]+"\n")
          pte.write(contents[i+1]+"\n")
          te_count += 1
    print("positive training examples : " + str(tr_count))
    print("positive validation examples : " + str(va_count))
    print("positive test examples : " + str(te_count))
    # neg eval
    ntr_count, nva_count, nte_count = 0, 0, 0

    for i in range(0, len(eval4contents), 2):
      if "N" not in eval4contents[i+1].upper():
        if eval4contents[i].split(">chr")[1].split(":")[0] in train_chroms and ntr_count < tr_count*ratio:
          ntr.write(eval4contents[i]+"\n")
          ntr.write(eval4contents[i+1]+"\n")
          ntr_count += 1
        elif eval4contents[i].split(">chr")[1].split(":")[0] in val_chroms and nva_count < va_count*ratio:
          nva.write(eval4contents[i]+"\n")
          nva.write(eval4contents[i+1]+"\n")
          nva_count += 1
        elif eval4contents[i].split(">chr")[1].split(":")[0] in test_chroms and nte_count < te_count*ratio:
          nte.write(eval4contents[i]+"\n")
          nte.write(eval4contents[i+1]+"\n")
          nte_count += 1
    print(ntr_count, nva_count, nte_count)
    # negative set

    for i in range(0, len(ncortex), 2):
      if "N" not in ncortex[i+1].upper():
        if ncortex[i].split(">chr")[1].split(":")[0] in train_chroms and ntr_count < tr_count*ratio:
          ntr.write(ncortex[i]+"\n")
          ntr.write(ncortex[i+1]+"\n")
          ntr_count += 1
        elif ncortex[i].split(">chr")[1].split(":")[0] in val_chroms and nva_count < va_count*ratio:
          nva.write(ncortex[i]+"\n")
          nva.write(ncortex[i+1]+"\n")
          nva_count += 1
        elif ncortex[i].split(">chr")[1].split(":")[0] in test_chroms and nte_count < te_count*ratio:
          nte.write(ncortex[i]+"\n")
          nte.write(ncortex[i+1]+"\n")
          nte_count += 1

    print(ntr_count, nva_count, nte_count)
  
    for i in range(0, len(ncontents), 2):
      if "N" not in ncontents[i+1].upper():
        if ncontents[i].split(">chr")[1].split(":")[0] in train_chroms and ntr_count < tr_count*ratio:
          ntr.write(ncontents[i]+"\n")
          ntr.write(ncontents[i+1]+"\n")
          ntr_count += 1
        elif ncontents[i].split(">chr")[1].split(":")[0] in val_chroms and nva_count < va_count*ratio:
          nva.write(ncontents[i]+"\n")
          nva.write(ncontents[i+1]+"\n")
          nva_count += 1
        elif ncontents[i].split(">chr")[1].split(":")[0] in test_chroms and nte_count < te_count*ratio:
          nte.write(ncontents[i]+"\n")
          nte.write(ncontents[i+1]+"\n")
          nte_count += 1
    print(ntr_count, nva_count, nte_count)
  ptr.close()
  pva.close()
  pte.close()
  ntr.close()
  nva.close()
  nte.close()


def make_mouse_splits(model_path, ratio):
  print("splitting mouse data...")
  with open(mouse_pos, "r") as f:
    pcontents = f.read().splitlines()
  f.close()
  with open(mouse_neg, "r") as f:
    ncontents = f.read().splitlines()
  f.close()
  with open(mouse_neg_cortex, "r") as f:
    ncortex = f.read().splitlines()
  f.close()
  with open(mouse_human_spec, "r") as f:
    nhuman = f.read().splitlines()
  f.close()

  # shuffle neg seqs
  negLocs = list()
  negSeqs = list()
  for i in range(0, len(ncontents), 2):
    if "N" not in ncontents[i+1].upper():
      negLocs.append(ncontents[i])
      negSeqs.append(ncontents[i+1])

  idxs = np.arange(len(negLocs))
  np.random.shuffle(idxs)

  with open(model_path+"/pos_TRAINING.fa", "w") as ptr, open(model_path+"/neg_TRAINING.fa", "w") as ntr,   open(model_path+"/pos_VALIDATION.fa", "w") as pva, open(model_path+"/neg_VALIDATION.fa", "w") as nva, open(model_path+"/pos_TEST.fa", "w") as pte, open(model_path+"/neg_TEST.fa", "w") as nte:
    # positive sequences
    tr_count, va_count, te_count = 0, 0, 0
    for i in range(0, len(pcontents), 2):
      if "N" not in pcontents[i+1].upper():
        if pcontents[i].split(":")[0].split(">chr")[1] in train_chroms:
          ptr.write(pcontents[i]+"\n")
          ptr.write(pcontents[i+1]+"\n")
          tr_count += 1
        elif pcontents[i].split(":")[0].split(">chr")[1] in val_chroms:
          pva.write(pcontents[i]+"\n")
          pva.write(pcontents[i+1]+"\n")
          va_count += 1
        elif pcontents[i].split(":")[0].split(">chr")[1] in test_chroms:
          pte.write(pcontents[i]+"\n")
          pte.write(pcontents[i+1]+"\n")
          te_count += 1
    print("positive training examples : " + str(tr_count))
    print("positive validation examples : " + str(va_count))
    print("positive test examples : " + str(te_count))

    # negative sequences
    ntr_count, nva_count, nte_count = 0, 0, 0
    for i in range(0, len(ncortex), 2):
      if "N" not in ncortex[i+1].upper():
        if ncortex[i].split(">chr")[1].split(":")[0] in train_chroms and ntr_count < tr_count*ratio:
          ntr.write(ncortex[i]+"\n")
          ntr.write(ncortex[i+1]+"\n")
          ntr_count += 1
        elif ncortex[i].split(">chr")[1].split(":")[0] in val_chroms and nva_count < va_count*ratio:
          nva.write(ncortex[i]+"\n")
          nva.write(ncortex[i+1]+"\n")
          nva_count += 1
        elif ncortex[i].split(">chr")[1].split(":")[0] in test_chroms and nte_count < te_count*ratio:
          nte.write(ncortex[i]+"\n")
          nte.write(ncortex[i+1]+"\n")
          nte_count += 1
        # negative sequences
    ntr_count, nva_count, nte_count = 0, 0, 0
    for i in range(0, len(nhuman), 2):
      if "N" not in "".join(nhuman[i+1]).upper():
        if nhuman[i].split("_")[0].split(">chr")[1] in train_chroms and ntr_count < tr_count*ratio:
          ntr.write(nhuman[i]+"\n")
          ntr.write("".join(nhuman[i+1])+"\n")
          ntr_count += 1
        elif nhuman[i].split("_")[0].split(">chr")[1] in val_chroms and nva_count < va_count*ratio:
          nva.write(nhuman[i]+"\n")
          nva.write("".join(nhuman[i+1])+"\n")
          nva_count += 1
        elif nhuman[i].split("_")[0].split(">chr")[1] in test_chroms and nte_count < te_count*ratio:
          nte.write(nhuman[i]+"\n")
          nte.write("".join(nhuman[i+1])+"\n")
          nte_count += 1
    for i in idxs:
      if "N" not in negSeqs[i].upper():
        if negLocs[i].split(":")[0].split(">chr")[1] in train_chroms and ntr_count < tr_count*ratio:
          ntr.write(negLocs[i]+"\n")
          ntr.write(negSeqs[i]+"\n")
          ntr_count += 1
        elif negLocs[i].split(":")[0].split(">chr")[1] in val_chroms and nva_count < va_count*ratio:
          nva.write(negLocs[i]+"\n")
          nva.write(negSeqs[i]+"\n")
          nva_count += 1
        elif negLocs[i].split(":")[0].split(">chr")[1] in test_chroms and nte_count < te_count*ratio:
          nte.write(negLocs[i]+"\n")
          nte.write(negSeqs[i]+"\n")
          nte_count += 1
    print("negative training examples : " + str(ntr_count))
    print("negative validation examples : " + str(nva_count))
    print("negative test examples : " + str(nte_count))


  ptr.close()
  pva.close()
  pte.close()
  ntr.close()
  nva.close()
  nte.close()


def main():
  path = "/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/"
  parser = argparse.ArgumentParser()
  parser.add_argument("-m", "--mouse", type=str, help="mouse model name")
  parser.add_argument("-hu", "--human", type=str, help="human model name")
  parser.add_argument("-r", "--ratio", type=float, help="ratio of negative to positive examples (maximum value is the nfold value used to generate negatives)")
  options, args = parser.parse_known_args()
  if (len(sys.argv)==1):
      parser.print_help(sys.stderr)
      sys.exit(1)
  elif (options.mouse is None and options.human is None):
    parser.print_help(sys.stderr)
    sys.exit(1)
  elif (options.ratio is None or options.ratio < 0):
    parser.print_help(sys.stderr)
    sys.exit(1)
  else:
    # make model directory
    if options.mouse is not None:
      os.makedirs(path+options.mouse)
      make_mouse_splits(path+options.mouse, options.ratio)
    if options.human is not None:
      os.makedirs(path+options.human)
      make_human_splits(halper_file, path+options.human, options.ratio)


if __name__ == "__main__":
  main()
