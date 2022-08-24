#!/bin/bash

# proces human retina samples 1 to 7
sbatch -p pfen_bigmem -J retina --export=ALL --mem 3G -t 4-0 --wrap "caper run ~/atac-seq-pipeline/atac_fix.wdl -i /projects/pfenninggroup/machineLearningForComputationalBiology/retina/scripts/zoonomia/human_ret.json"

# process human retina samples 2, 3, 4, 5, 6, 7
#sbatch -p pfen1 -J human_retina --export=ALL --mem 3G -t 4-0 --wrap "caper run ~/atac-seq-pipeline/atac_fix.wdl -i /projects/pfenninggroup/machineLearningForComputationalBiology/retina/scripts/zoonomia/human_ret.json"

# process mouse retina samples 1, 2, 3
#sbatch -p pfen1 -J mouse_retina --export=ALL --mem 3G -t 4-0 --wrap "caper run ~/atac-seq-pipeline/atac_fix.wdl -i /projects/pfenninggroup/machineLearningForComputationalBiology/retina/scripts/zoonomia/mouse_ret_GSE153817.json"

# process mouse retina samples 1, 2, 3
#sbatch -p pfen1 -J mouse_retina --export=ALL --mem 3G -t 4-0 --wrap "caper run ~/atac-seq-pipeline/atac_fix.wdl -i /projects/pfenninggroup/machineLearningForComputationalBiology/retina/scripts/zoonomia/mouse_ret_GSE146897.json"

