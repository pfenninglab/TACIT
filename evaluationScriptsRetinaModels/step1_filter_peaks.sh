#!/bin/bash

# human retina atac-seq

HUMAN_PATH="/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/human/GSE137311/atac/6a26719a-12e5-4901-bd59-f0c4c5cd9495/call-reproducibility_idr/execution/idr.optimal_peak.narrowPeak.gz"
echo "human IDR optimal peaks : "
zcat ${HUMAN_PATH} | wc -l

# # remove optimal IDR peaks within 20 kb of a TSS
NO_TSS_PATH="/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/human/GSE137311/hg38_ret_noTSS.bed"
bedtools window -v -w 20000 -a $HUMAN_PATH -b ~/resources/hg38_gencode_tss_unique.bed.gz > $NO_TSS_PATH
echo "human no TSS : "
cat $NO_TSS_PATH | wc -l

# # subtract blacklist regions
# # source : https://github.com/Boyle-Lab/Blacklist/tree/master/lists
echo "human no TSS, no blacklist : "
bedtools window -v -a $NO_TSS_PATH -b ~/resources/hg38-blacklist.v2.bed.gz | wc -l

# # remove protein coding exons, any peaks longer than 1 kb
FILTERED_PATH="/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/human/GSE137311/hg38_ret_noTSS_filtered.bed"
bedtools intersect -v -a $NO_TSS_PATH -b /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanGenome/gencode.v27.annotation.protExon_geneNames.bed | bedtools window -v -a stdin -b ~/resources/hg38-blacklist.v2.bed.gz | awk '{if ($3-$2 <= 1000) {print $0}}' > $FILTERED_PATH
echo "human no TSS, no blacklist, no protein coding exons, no long peak: "
cat $FILTERED_PATH | wc -l

# mouse retina atac-seq
#echo ""
#MOUSE_PATH="/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/mouse/GSE153817_WTMouse_ret_ATAC/atac/05006732-b3c9-45fb-bd43-8af86e6e4765/call-reproducibility_idr/execution/idr.optimal_peak.narrowPeak.gz"
#MOUSE_PATH="/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/mouse/GSE146897_WTMouse_ret_ATAC/atac/9c087798-a245-403d-a2d0-330d2c2dd7cf/call-reproducibility_idr/execution/idr.optimal_peak.narrowPeak.gz"
#echo "mouse IDR optimal peaks : "
#zcat ${MOUSE_PATH} | wc -l

# remove optimal IDR peaks within 20 kb of a TSS
#NO_TSS_PATH="/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/mouse/GSE146897_WTMouse_ret_ATAC//mm10_ret_noTSS.bed"
#bedtools window -v -w 20000 -a $MOUSE_PATH -b ~/resources/mm10TSS_RefSeq.bed > $NO_TSS_PATH
#echo "mouse no TSS : "
#cat $NO_TSS_PATH | wc -l

# subtract blacklist regions
# source : https://github.com/Boyle-Lab/Blacklist/tree/master/lists
#echo "mouse no TSS, no black list: "
#bedtools window -v -a $NO_TSS_PATH -b ~/resources/mm10-blacklist.v2.bed.gz | wc -l

# remove any peaks longer than 1 kb
#FILTERED_PATH="/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/mouse/GSE146897_WTMouse_ret_ATAC//mm10_ret_noTSS_filtered.bed"
#bedtools window -v -a $NO_TSS_PATH -b ~/resources/mm10-blacklist.v2.bed.gz | bedtools window -v -a stdin -b ~/resources/mm10_protein_coding_exons.bed | awk '{if ($3-$2 <= 1000) {print $0}}' > $FILTERED_PATH
#echo "mouse no TSS, no blacklist, no protein coding exons, no long peaks : "
#cat $FILTERED_PATH | wc -l
