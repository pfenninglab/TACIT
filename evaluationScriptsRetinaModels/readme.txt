
File Descriptions:

step0_process_atacseq.sh : A shell script that launches the ATAC-seq pipeline to call peaks on raw, bulk human and mouse retina ATAC-seq data

step1_filter_peaks.sh : A shell script that performs the following filters on the retina optimal IDR peaks: removes peaks within 20 kB of a transcription start site, removes peaks overlapping ENCODE blacklist v2 regions, removes peaks overlapping protein coding exons, and removes peaks longer than 1 kB in length.

step2_label_peaks.sh : A shell script that labels filtered peak and summit BED files with unique numbers as a preprocessing step for ortholog mapping.

step3_run_peaks_hal.sb : A sbatch script that submits a job to map peaks from one species references genome to another.

step4_run_summits_hal.sb : A sbatch script that submits a job to map peak summits from one species reference genome to another.

step5_run_halper.sh : A shell script that runs HALPER on the mapped peaks and summits.

step6a_generate_neg_seqs_mm10.sb : A sbatch script that generates 10X G-C matched negative sequences for mouse retina filtered open chromatin sequences and outputs the sequences in a FASTA file.

step6b_generate_neg_seqs_hg38.sb : A sbatch script that generates 10X G-C matches negative sequences for human retina filtered open chromatin sequences and outputs the sequences in a FASTA file.

step7_split_trainvaltest.py : A Python 3 script that partitions the retina filtered open chromatin sequences by chromosome into training, validation, and test sets as input datasets for machine learning models.

step9c_keras_cnn.py : A Python 3 script that trains a convolutional neural network to classify DNA sequences underlying open chromatin by their cell type membership.

step10_roc_evaluations.py : A Python 3 script that produces visualizations of the convolutional neural network performance on data that the model was not trained on.

step11c_predict_glires.py : A Python 3 script that uses the convolutional neural network to annotate orthologs of retinal open chromatin sequences,

step12a_deepshap.py : A Python 3 script that uses the DeepShAP algorithm to score retinal open chromatin sequences for their per-base importance of contribution towards model classification.

step12b_modisco.py : A Python 3 script that uses the DeepShAP scores to generate position weight matrices representing transcription factor motifs that were learned by the convolutional neural network.

step14b_report_metrics.py : A Python 3 script that ouputs all reported test set metrics.

tfmodisco_visualize.ipynb : A iPython notebook that shows visualizations corresponding to the DeepShAP and TF-MoDISco analyses.

keras2.yml : An anaconda Python 3 environment containing all the necessary packages to run the Python 3 scripts in this folder.
