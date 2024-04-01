# TACIT Enhancer-Phenotype Association
Enhancer-Phenotype Association from Tissue-Aware Conservation Inference Toolkit.

We used this with open chromatin regions (OCRs) / predictions of open chromatin status (1), but the pipeline is also amenable to predictions of other proxies for enhancer activity.  This pipeline begins with the step of TACIT for associating enhancer activity predictions with phenotypes (step 4 of Figure 1 from Kaplow*, Lawler*, Schaffer*, _et al_., Science, 2023 (1)) and continues through the end of the pipeline; an explanation of how to implement steps 1-3 can be found near the end of this README, and details can be found in this document: https://docs.google.com/document/d/1_fVGWQ6UYhM4oBszjCu04upOu9MzQhwEzHoLJ1n2iYQ/edit?usp=sharing.

This pipeline has been tested on Centos versions 7 and 8. However, we expect it to work on other distributions of Linux, MacOS, and Windows.


## Recommendations for computing enhancer-phenotype associations
1. Run `ocr_phylo(g)lm.r` (described below).  If there are no OCRs with p-value < 0.05, then predicted open chromatin in the tissue or cell type of interested is not associated with the phenotype, so stop.
2. Run `permulationList.py` (described below) on the output from `ocr_phylo(g)lm.r` with --threshold 1 and --permulation 999.
3. Run `ocr_phylo(g)lm_conditional.r` (described below) on the output from `permulationList.py`.
4. Run `compute_perm_pvals_conditional.py` (described below) on the output from `ocr_phylo(g)lm_conditional.r`.
5. Run `permulationList.py` on the output from `compute_perm_pvals_conditional.py` with threshold 0.05 and --permulation 9000.
6. If OCRs remain, repeat steps 3-5, using threshold 0.005 and --permulation 90000 in step 5.
7. If OCRs remain, repeat steps 3-5, using threshold 0.0005 and --permulation 900000 in step 5.
8. If OCRs remain, repeat steps 3-4.
9. Run `bhCorrection.R` (described below) on the output from `compute_perm_pvals_conditional.py`.


## Programs for OCR-phenotype association with phylo(g)lm (2) and phylogenetic permulations (3):
`ocr_phylolm.r` and `ocr_phyloglm.r`: Perform OCR-phenotype associations for continuous and binary traits, using phylolm and phyloglm, respectivly. Can generate unpermulated results as well as a fixed number of permulations for each OCR. See the header comment of this program for usage details, especially for external parallelization.

_Input arguments (must be provided in order)_:
1. tree file (Newick file)
2. open chromatin predictions matrix file (text file, rows are OCRs, columns are species no header, first column is OCR names, tab-seperated, missing data is noted as -1; examples can be found at http://daphne.compbio.cs.cmu.edu/files/ikaplow/TACITSupplement/PredictionMatrices/)
3. file listing species corresponding to predictions matrix columns (text file, one per line, in the format Genus species)
4. phenotype file (csv file with column name headers, species names in a column called "species binomial" in the format genus_species for `ocr_phylolm.r` and  "Species Name" in the format Genus species for `ocr_phyloglm.r`)
5. output file name "template," which includes a directory followed by a csv file name; output files will be the file name without the csv followed by the output file number and the random seed (example: if output file "template" given is /path/to/foo.csv, output will be in (I,S as below) /path/to/foo_rI_S.csv)
6. I: output file number and inital line in predictions matrix
7. J: step size in predictions matrix (see below),
8. K: number of permulations per OCR (0 for true data / no permulations, make > 0 only if not using additional rejection sampler for accepting only permulations that preserve the direction from phylolm/phyloglm as described below; will apply phylolm to K permulations for OCRs on lines I, I+J, I+2J, ... until end of matrix is reached)
9. S: random seed to use
10. [ONLY FOR `ocr_phyloglm.r`] path to directory containing `fast_bin_perm.r`
11. name of column in phenotype file with phenotype
12. (optional) name of column with phenotype that should be treated as a covariate (can be more than one)

Example (submitted to slurm cluster in sbatch script): Rscript ocr_phylolm.r Zoonomia_ChrX_lessGC40_241species_30Consensus.tree cortex_ocrs_filtered_named.txt species.txt Zoonomia_phenotypes_12-14-21.csv BrainSizeResidMotorCortexPhylolmOut/motorCortex_brainSizeResid_results.csv ${SLURM_ARRAY_TASK_ID} 1 0 1 Brain.resid
                  
`ocr_phylolm_conditional.r` and `ocr_phyloglm_conditional.r`: Conditional-test versions of the above scripts. Both have an additrional rejection sampler that only accepts permulations that associate with activity predictions in the same direction as the true phenotype. Instead of a fixed number of trials, both expect as input an CSV containing the original (unpermulated) coefficient and the number of trials to perform for each OCR. The specification for this input file matches the output of the conditional scripts for computing p-values, described below.

_Input arguments (must be provided in order)_:
1. tree file (Newick file)
2. open chromatin predictions matrix file (text file, rows are OCRs, columns are species no header, first column is OCR names, tab-seperated, missing data is noted as -1; examples can be found at http://daphne.compbio.cs.cmu.edu/files/ikaplow/TACITSupplement/PredictionMatrices/)
3. file listing species corresponding to predictions matrix columns (text file, one per line, in the format Genus species)
4. phenotype file (csv file with column name headers, species names in a column called "species binomial" in the format genus_species for `ocr_phylolm.r` and  "Species Name" in the format Genus species for `ocr_phyloglm.r`)
5. output file name "template," which includes a directory followed by a csv file name; output files will be the file name without the csv followed by the output file number and the random seed (example: if output file "template" given is /path/to/foo.csv, output will be in (I,S as below) /path/to/foo_rI_S.csv)
6. I: output file number and inital line in predictions matrix
7. J: step size in predictions matrix
8. file with columns OCR that contains OCR names and Missing_Trials that specifies how many permulations (K) to do for each OCR (csv file, can make by running `permulationList.py` on output from `compute_perm_pvals_conditional.py` below)
9. S: random seed to use
10. [ONLY FOR `ocr_phyloglm_conditional.r`] path to directory containing `fast_bin_perm.r`
11. name of column in phenotype file with phenotype
12. (optional) name of column with phenotype that should be treated as a covariate (can be more than one)

Example (submitted to slurm cluster in sbatch script): Rscript ocr_phylolm_conditional.r Zoonomia_ChrX_lessGC40_241species_30Consensus.tree cortex_ocrs_filtered_named.txt species.txt Zoonomia_phenotypes_12-14-21.csv BrainSizeResidMotorCortexPhylolmPermulationsOut/motorCortex_brainSizeResid_results_round1.csv ${SLURM_ARRAY_TASK_ID} 1000 motorCortex_brainSizeResid_results_modified_peakList.csv 1 Brain.resid

`fast_bin_perm.r`: Our fast implementation of binary permulations (1). It contains the function `fastSimBinPhenoVec`, which is called by `ocr_phyloglm.r` and can also be used to generate permulations of binary traits for other statistical tests. The function `makeLeafMap` precomputes the required data structure. The other functions are called only internally.


## Utilities:
`compute_perm_pvals.py`: Script for computing permulation p-values using the output from `ocr_phylo(g)lm.r`. It takes an input a file containing the unpermulated results (output of `ocr_phylo(g)lm.r`) and a directory containing permulated results. See the header comment of this script for usage details.  We recommend using this method only when compute time is limited and the method below when sufficient compute time is available.

_Input arguments (must be provided in order)_:
1. the result file of `ocr_phylo(g)lm.r` run once without permulations for every OCR (csv file)
2. the number of directories containing results files from permulations
3. directories containing results file(s) of ocr_phylo(g)lm run with permulations on a subset of the OCRs in the first file (the number of directories provided should equal the previous argument)
4. an output csv file, which will be a copy of the input file with additional columns giving the permulations p-value and number of trials performed for each OCR

Example: python compute_perm_pvals.py BrainSizeResidMotorCortexPhylolmOut/motorCortex_brainSizeResid_results_r1_s1.csv 1 BrainSizeResidMotorCortexPhylolmPermulationsOut motorCortex_brainSizeResid_results_perm1k.csv

`compute_perm_pvals_conditional.py` is an alternate version of compute_perm_pvals.py that, instead of using a pure one-sided test, rejects permulations in the opposite direction. This is the conditional test described in (4) with A=0 for all OCRs. We recommend using this method instead of compute_perm_pvals.py for traits that can be easily permulated, as it is more principled but slower due to its rejection of permulations in the opposite direction. It otherwise behaves the same as compute_perm_pvals.py, except that this script allows for multiple input directories (preceded by a count). Additionally, the output contains a `Missing_Trials` column, which is a guess of how many permulations need to be replaced due to rejection. When used as input to `ocr_phylo(g)lm_conditional.r`, this column also specifies the number of trials to perform foreach OCR, and manual editing may be needed to get the desired number of permulations.

_Input arguments (must be provided in order)_:
1. the result file of `ocr_phylo(g)lm.r` run once without permulations for every OCR (csv file)
2. the number of directories containing results files from permulations
3. directories containing results file(s) of ocr_phylo(g)lm run with permulations on a subset of the OCRs in the first file (the number of directories provided should equal the previous argument)
4. an output csv file, which will be a copy of the input file with additional columns giving the permulations p-value and number of trials performed for each OCR

Example: python compute_perm_pvals_conditional.py BrainSizeResidMotorCortexPhylolmOut/motorCortex_brainSizeResid_results_r1_s1.csv 1 BrainSizeResidMotorCortexPhylolmPermulationsOut motorCortex_brainSizeResid_results_perm1k.csv

`compute_perm_pvals_conditional_extended.py` is a variant of the above script that offers some additional functionality. First, if an OCR occurs in the permulated results but not the original results file, it is skipped rather than raising an error. Second, it takes additional positional arguments specifying a maximum number of permulations and an additional output file; if more permulations than the maximum number are found, they are written to the additional output file for possible use in a later round. Third, the number of "better" trials (those in the numerator of the p-value) is explicitly reported as the 7th column of the output. And, if the "original" results file contains at least seven columns, it is assumed to be the output of this script and the number of "better" and total trials are initialized to those in the file, rather than 1. Fourth, if the unpermulated coefficient of an OCR is exactly 0, permulations whose coefficient is also exactly 0 will be accepted; this will not affect p-values in this case (which will always be 1) but will affect the reported number of trials for those OCRs.

`permulationList.py`: Prepare output from phylogenetic permulations for an additional round of phylogenetic permulations.  This updates the "Missing_Trials" column to add additional missing trials for the subset of rows with sufficiently low p-values to warrent additional permulations.

_Input arguments_:
1. --input: csv file with Exp_Pvalue column and a column with peak name
2. --threshold:  threshold used to filter peaks by p-values (remove peaks that are not close to significant for future rounds of permulations, use 1 if no permulations have been run)
3. --permulation: value to put in Missing_Trials column (number of additional permulations to do)
4. --matrix: open chromatin predictions matrix file (text file, rows are OCRs, columns are species no header, first column is OCR names, tab-seperated, missing data is noted as -1; examples can be found at http://daphne.compbio.cs.cmu.edu/files/ikaplow/TACITSupplement/PredictionMatrices/)
5. --colName: name of column name with p-values used for filtering (Pvalue if no permulations have been run, Exp_Pvalue if permulations have been run)
6. --output: output file name prefix, where one output file will be file with columns OCR that contains OCR names and Missing_Trials that specifies how many permulations to do for each OCR (csv file with suffix peakList.csv) and the other will be a the open chromatin matrix file with the subset of lines with p-value lower than the threshold

Example: python permulationList.py --input motorCortex_brainSizeResid_results_perm1k.csv --threshold 0.05 --permulation 9000 --matrix cortex_ocrs_filtered_named.txt --colName Exp_Pvalue --output motorCortex_brainSizeResid_results_perm1k_pLessThan0.05

`bhCorrection.R`: Run the Benjamini-Hochberg Procedure (5) on the p-values from phylogenetic permulations and create an output file that has the adjusted p-values, sorted from lowest to highest.

_Input arguments (must be provided in order)_:
1. input file with OCR names as well as the p-values from phylo(g)lm, coefficients from pylo(g)lm, permulations p-values, number of trials performed for each OCR (csv file, obtained from `compute_perm_pvals.py` or `compute_perm_pvals_conditional.py`)
2. output file name, where output file will be input file with additional column with the Benjamini-Hochberg adjusted p-values

Example: Rscript /home/ikaplow/RegulatoryElementEvolutionProject/src/TACIT/bhCorrection.R motorCortex_brainSizeResid_results_perm1m.csv motorCortex_brainSizeResid_results_perm1m_adjp.csv


## Dependencies:
### For R programs:
R >= 4.0

ape package >= 5.4.1 (6)

phylolm package >= 2.6.2 (2)

geiger package >= 2.0.7 (7)

Rcpp package >= 1.0.6 (8)

stringr package >= 1.4.0


### For Python programs:
Python 3

pandas 0.2.24-1.3.5


## Scripts for generating figures in Kaplow, Lawler, Schaffer, _et al_. (1):

evaluationScriptsRERconverge: scripts for preparing data for and running RERconverge

evaluationScriptsCortexLiverModels: scripts for evaluating and making predictions with motor cortex and liver models

evaluationScriptsPVModels: scripts for evaluating and making predictions with PV models

evaluationScriptsRetinaModels: scripts for evaluating and making predictions with retina models

evaluationScriptsEnhancerPhenotypeAssociations: scripts for evaluating and analyzing open chromatin region-phenotype associations

evaluationScriptsHiC: scripts for evaluating long-range interactions between phenotype-associated open chromatin regions and phenotype-relevant genes

evaluationScriptsVocalLearning: scripts for evaluating associations between open chromatin regions and vocal learning and integrating results with gene-vocal learning associations

Dependencies specific to those scripts are listed in the individual directories.  Models can be found at http://daphne.compbio.cs.cmu.edu/files/ikaplow/TACITSupplement/.


## Instructions for obtaining predictions from models trained in Kaplow, Lawler, Schaffer, _et al_. (1):
1. Go to the UCSC Genome Browser (genome.ucsc.edu).
2. Under "My Data," select "Track Hubs."
3. This will take you to the Track Data Hubs pages, which has three tabs.  Select "Connected Hubs."
4. Paste in https://cgl.gi.ucsc.edu/data/cactus/241-mammalian-2020v2-hub/hub.txt, and click "Add Hub."
5. For obtaining the predictions in the future, go to https://genome.ucsc.edu/cgi-bin/hgGateway?genome=Homo_sapiens&hubUrl=https://cgl.gi.ucsc.edu/data/cactus/241-mammalian-2020v2-hub/hub.txt, where the genome= parameter can be replaced with other species.


## Implementing steps 1-3 of figure 1 from Kaplow*, Lawler*, Schaffer*, _et al_., Science, 2023 (1):

Step 1: For processing bulk open chromatin data, we used the ENCODE open chromatin pipeline (https://github.com/ENCODE-DCC/atac-seq-pipeline).  For processing single-cell ATAC-seq data (step 1), we used Archr (https://www.archrproject.com/).  We removed potential promoters (OCRs within 20kb of a protein-coding transcription start site), potential super-enhancers (OCRs longer than 1kb), and OCRs overlapping protein-coding exons.

Step 2: For training machine learning models to use DNA sequence to predict open chromatin (step 2), we used keras (https://keras.io/).  Our machine learning models can be found at http://daphne.compbio.cs.cmu.edu/files/ikaplow/TACITSupplement/CNNs/.  The models for predicting brain open chromatin and the models for predicting liver open chromatin trained using only mouse sequences can be found at https://github.com/pfenninglab/OCROrthologPrediction/tree/main/models (brain model architechture: brainEnhancer_humanMouseMacaqueRat_euarchontaglireEnhLooseOrthNeg_500bp_conv5_architecture.json, brain model weights: brainEnhancer_humanMouseMacaqueRat_euarchontaglireEnhLooseOrthNeg_500bp_conv5.hdf5, liver model architecture: liverEnhancer_euarchontaglireEnhLooseOrthNeg_500bp_conv5_architecture.json, liver model weights: liverEnhancer_euarchontaglireEnhLooseOrthNeg_500bp_conv5.hdf5).  For comparing results to DeepSEA Beluga, we input our sequences into https://humanbase.net/deepsea/ with the Beluga model option.

Step 3: For mapping OCR orthologs across species, we used halLiftover (https://github.com/ComparativeGenomicsToolkit/hal) followed by HALPER (https://github.com/pfenninglab/halLiftover-postprocessing).  For making predictions at orthologs we used `predictNewSequencesNoEvaluation.py` from https://github.com/pfenninglab/OCROrthologPrediction.


## Contact
Irene Kaplow (ikaplow@cs.cmu.edu)

Daniel Schaffer (dschaffe@andrew.cmu.edu)

Alyssa Lawler (alawler@andrew.cmu.edu)

Heather Sestili (hharper@cmu.edu)

Chaitanya Srinivasan (csriniv1@andrew.cmu.edu)

Tianyu Lu (tianyul3@andrew.cmu.edu)

Andreas Pfenning (apfenning@cmu.edu)


## References
1. I. M. Kaplow, A. J. Lawler, D. E. Schaffer, C. Srinivasan, H. H. Sestili, M. E. Wirthlin, B. N. Phan, K. Prasad, A. R. Brown, X. Zhang, K. Foley, D. P. Genereux, Zoonomia Consortium, E. K. Karlsson, K. Lindblad-Toh, W. K. Meyer, A. R. Pfenning, Relating enhancer genetic variation across mammals to complex phenotypes using machine learning. *Science*. **380**, eabm7993 (2023).
2. L. S. T. Ho, C. Ané, A linear-time algorithm for Gaussian and non-Gaussian trait evolution models. *Syst. Biol*. **63**, 397–408 (2014).
3. E. Saputra, A. Kowalczyk, L. Cusick, N. Clark, M. Chikina, Phylogenetic Permulations: A Statistically Rigorous Approach to Measure Confidence in Associations in a Phylogenetic Context. *Mol. Biol. Evol*. **38**, 3004–3021 (2021).
4. E. Kulinskaya, On two-sided p-values for non-symmetric distributions. *arXiv*. https://arxiv.org/abs/0810.2124 (2008).
5. Y. Benjamini and Y. Hochberg, Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. *Journal of the Royal Statistical Society*. **57**, 289-300 (1995).
6. E. Paradis, K. Schliep, ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. *Bioinformatics*. **35**, 526–528 (2019).
7. M. W. Pennell, J. M. Eastman, G. J. Slater, J. W. Brown, J. C. Uyeda, R. G. FitzJohn, M. E. Alfaro, L. J. Harmon, geiger v2.0: an expanded suite of methods for fitting macroevolutionary models to phylogenetic trees. *Bioinformatics*. **30** (2014), pp. 2216–2218.
8. D. Eddelbuettel, R. François, J. Allaire, K. Ushey, Q. Kou, N. Russel, J. Chambers, D. Bates, Rcpp: Seamless R and C++ integration. *J. Stat. Softw*. **40**, 1–18 (2011).
