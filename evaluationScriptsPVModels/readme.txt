step1_make_BICCN_mouse_NEMO_ArchRProject.r: makes an ArchR project for the mouse motor cortex snATAC-seq data from Li et al., Nature, 2021 (https://www.nature.com/articles/s41586-021-03604-1)

step2_plot_BICCN_mouse_NEMO_ArchRProject.r: plots heatmaps for mouse motor cortex snATAC-seq data processed by ArchR

step3_extractPeaks_BICCN_mouse_NEMO_ArchRProject.r: gets reproducible peaks in neurons for mouse motor cortex snATAC-seq data processed by ArchR

step1_make_BICCN_huMOp_ArchRProject.r: makes an ArchR project for the human motor cortex snATAC-seq data from SNARE-seq2 from Bakken et al., Nature, 2021 (https://www.nature.com/articles/s41586-021-03465-8)

step2_plot_BICCN_huMOp_ArchRProject.r: plots heatmaps for human motor cortex snATAC-seq data from SNARE-seq2 processed by ArchR

step3_extractPeaks_BICCN_mouse_NEMO_ArchRProject.r: gets reproducible peaks in neurons for human motor cortex snATAC-seq data from SNARE-seq2 processed by ArchR

printCommandsFor_HALPER2bed.py: creates a shell script for converting mapped peak regions to bed files of 500bp summit centered orthologs for each species

printCommandsFor_getfasta.py: creates a shell script for extracting fasta sequences for the bed regions for each species from that species' genome 

generatePredictions_PVcnn.py: scores fasta sequences through the PV CNN for each species and outputs the scores to a text file

PVmodelEvaluations.py: gets performance metrics for evaluation criteria described in the paper
