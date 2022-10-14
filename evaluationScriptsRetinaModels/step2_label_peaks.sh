#!/bin/bash

# narrowPeak filtered peaks 
HUMAN_PEAKS="/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/human/GSE137311/hg38_ret_noTSS_filtered.bed"
MOUSE_PEAKS="/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/mouse/GSE153817_WTMouse_ret_ATAC/mm10_ret_noTSS_filtered.bed"

# label peaks 
bedtools sort -i $HUMAN_PEAKS | awk 'OFS="\t" {print $1, $2, $3, "peak"NR-1}' > "${HUMAN_PEAKS::-4}_labeled.bed"
#bedtools sort -i $MOUSE_PEAKS | awk 'OFS="\t" {print $1, $2, $3, "peak"NR-1}' > "${MOUSE_PEAKS::-4}_labeled.bed"
# flank 200 bp
awk 'OFS="\t" {print $1, $2-200, $3+200, $4}' > "${HUMAN_PEAKS::-4}_labeled_flank200.bed"

# define summits with labels
bedtools sort -i $HUMAN_PEAKS | awk 'OFS="\t" {print $1, $2+$10, $2+$10+1, "peak"NR-1}' > "${HUMAN_PEAKS::-4}_summits.bed"
#bedtools sort -i $MOUSE_PEAKS | awk '{print $1, $2+$10, $2+$10+1, "peak"NR-1}' > "${MOUSE_PEAKS::-4}_summits.bed"