#!/bin/bash

###################### HUMAN #################################

# run halper on halLiftover outputs
FILE="/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/human/GSE137311/hg38_ret_noTSS_filtered_"

python ~/halLiftover-postprocessing/orthologFind.py -max_len 1000 -min_len 50 -protect_dist 5 -max_frac 5.0 -min_frac 0.05 -qFile "${FILE}labeled_500bp.bed" -tFile "${FILE}labeled_halLiftover_mm10.bed" -sFile "${FILE}summits_halLiftover_mm10.bed" -oFile "${FILE}halper_mm10.bed" -mult_keepone
OLD=$(< "${FILE}labeled_500bp.bed" wc -l )
NEW=$(< "${FILE}halper_mm10.bed" wc -l )
RESULT=$(echo "$NEW/$OLD" | bc -l)
echo "mappability : " $RESULT

# # retain only mappable peaks
awk -F "\t" 'NR==FNR{F1[$5];next}$4 in F1{print}' "${FILE}halper_mm10.bed" "${FILE}labeled.bed" > "${FILE}mappable_peaks.bed"

# # get mouse (mm10) sequences and human-mouse (hg38) mappable sequences for MEME-ChIP
# # sequences must be equal (500 bp) and summit centered

# # match mappable peaks to their summits and take 500 bp sequences
awk -F "\t" 'NR==FNR{F1[$4];next}$4 in F1{print}' "${FILE}mappable_peaks.bed" "${FILE}summits.bed" | awk 'OFS="\t" {print $1, $2-250, $2+250, $4}' > "${FILE}mappable_peaks_500bp.bed"
#awk 'OFS="\t" {print $1, $2-250, $2+250, $4}' /projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/mouse/GSE153817_WTMouse_ret_ATAC/mm10_ret_noTSS_filtered_summits.bed > /projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/mouse/GSE153817_WTMouse_ret_ATAC/mm10_ret_noTSS_filtered_500bp.bed

bedtools getfasta -fi /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanGenome/hg38.fa -bed "${FILE}mappable_peaks_500bp.bed" > "${FILE}mappable_peaks_500bp.fa"
#bedtools getfasta -fi ~/resources/mm10.fa -bed /projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/mouse/GSE153817_WTMouse_ret_ATAC/mm10_ret_noTSS_filtered_500bp.bed > /projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/mouse/GSE153817_WTMouse_ret_ATAC/mm10_ret_noTSS_filtered_500bp.fa


################################### MOUSE ###########################################s

#FILE="/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/mouse/GSE153817_WTMouse_ret_ATAC/mm10_ret_noTSS_filtered_"

#python ~/halLiftover-postprocessing/orthologFind.py -max_len 1000 -min_len 50 -protect_dist 5 -max_frac 5.0 -min_frac 0.05 -qFile "${FILE}labeled.bed" -tFile "${FILE}labeled_halLiftover_hg38.bed" -sFile "${FILE}summits_halLiftover_hg38.bed" -oFile "${FILE}halper_hg38.bed" -mult_keepone
#OLD=$(< "${FILE}labeled.bed" wc -l )
#NEW=$(< "${FILE}halper_hg38.bed" wc -l )
#RESULT=$(echo "$NEW/$OLD" | bc -l)
#echo "mappability : " $RESULT

# retain only mappable peaks
#awk -F "\t" 'NR==FNR{F1[$5];next}$4 in F1{print}' "${FILE}halper_hg38.bed" "${FILE}labeled.bed" > "${FILE}mappable_peaks.bed"

# match mappable peaks to their summits and take 500 bp sequences
#awk -F "\t" 'NR==FNR{F1[$4];next}$4 in F1{print}' "${FILE}mappable_peaks.bed" "${FILE}summits.bed" | awk 'OFS="\t" {print $1, $2-250, $2+250, $4}' > "${FILE}mappable_peaks_500bp.bed"
#bedtools getfasta -fi /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanGenome/hg38.fa -bed "${FILE}mappable_peaks_500bp.bed" > "${FILE}mappable_peaks_500bp.fa"
