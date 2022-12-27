Dependencies.txt lists dependencies for files in this directory beyond those listed in the main directory.

ocr_phylolm_figures.r contains code was used in the generation of figures 4, 5, 6C-D, S8, and S9. 
It includes some code from ocr_phylolm.r, for use in loading and generating data, as well as several blocks of
code below that were used to generate plots showing associations between OCR activity predictions and phenotypes.
The code included contains various hardcoded filenames and values that were changed to generate panels for different OCRs and phenotypes.
This script is not intended to be run as a whole; rather, different blocks of code were run to produce different figures.
It is provided to illustrate how figures were produced and is not intended for reuse. 

pdists.R is a small script used for plotting histograms of the p-value distributions for associations between (motor cortex and PV neuron) OCRs and (brain size and solitary living).

p_dist_collect.py is a Python 3 script for comparing p-value distributions for p-values associated with OCRs that are
near or not near certain genes. It can also be used for p-value-like metrics, such as the fraction of successful trials. It takes four positional command-line arguments:
1) A list of genes of interest, one per line
2) A space-separated file with each gene near each OCR, one pair per line, with the OCRs in column 4 and the genes in column 8
3) A space-separated file with the p-value associated with each OCR (one per line), with OCRs in column 1 and p-values in column 4.
4) A file to write output, which is a CSV where the first column contains the p-values of each OCR, and the second column contains
   T if the corresponding OCR is near a gene of interest and F otherwise.
The program writes output as described, and additionally reports the result of a Wilcoxon rank-sum test between the set of p-values
associated with OCRs near genes of interest and the set of p-values near other genes but not genes of interests. 
OCRs near no genes (i.e. that do not occur in column 4 of input file 2) are ignored.

p_dist_collect_gene.py is a Python 3 script for comparing p-value distributions for p-values associated with genes that are
near or not near OCRS. It takes three positional command-line arguments:
1) A space-separated file with each gene near each OCR, one pair per line, with the OCRs in column 4 and the genes in column 8
2) A space-separated file with the p-value associated with  each gene (one per line), with genes in column 1 and p-values in column 4.
3) A file to write output, which is a CSV where the first column contains the p-values of each gene, and the second column contains
   T if the corresponding gene is near an OCR and F otherwise.
The program writes output as described, and additionally reports the result of a Wilcoxon rank-sum test between the set of p-values
associated with genes near OCRs of interest and the set of p-values associated with all other genes that occur in input file 3. 
OCRs of interest are those that appear in column 4 of input file 1.

phylop_p_correlation is a Python 3 script for comparing p-value distributions for p-values associated with OCRs with average PhyloP conservation scores. It can also be used for p-value-like metrics, such as the fraction of successful trials. It takes three positional command-line arguments:
1) Mouse or Human coordinates of each OCR, as a BED file. 
2) A comma-separated file with the p-value associated with each OCR (one per line), with OCRs in column 1 and p-values in column 4.
3) [human | mouse], to specify the species to use for PhyloP scores. Paths for each species are hardcoded. 
The program reports the result of a Pearson correlation between the set of p-values associated with OCRs and the average PhyloP scores of those OCRs.

distance_p_correlation is a Python 3 script for comparing p-value distributions for p-values associated with OCRs with genemoic distances. It can also be used for p-value-like metrics, such as the fraction of successful trials. It takes three positional command-line arguments:
1) Mouse or Human coordinates of each OCR, as a BED file, with a distance (or other metric) as the last column
2) A comma-separated file with the p-value associated with each OCR (one per line), with OCRs in column 1 and p-values in column 4.
The program reports the result of a Pearson correlation between the set of p-values associated with OCRs and the distances associated with those OCRs.

makePGLSSbatch.py is a Python 3 script for creating a submission script for running ocr_phylolm.r or ocr_phyloglm.r on a cluster.

runPGLSPSCFilt.sh is a shell script that was used to associate binary phenotypes with OCR ortholog open chromatin predictions.

correctSolitaryResults.r is an R script that does Benjamini-Hochberg correction for solitary and group living associations.

correctVocalLearningResults.r is an R script that does Benjamini-Hochberg correction for vocal learning associations.

sim_perm_p.val_subset.R is an R script that computes t-test p-values for simulated data with different fractions of features from samples drawn from different distributions using the parametric t-test and the emperical permutation-based method used for TACIT that is akin to a one-tailed test.

sim_perm_p.val_newpVal.R is an R script that computes t-test p-values for simulated data with different fractions of features from samples drawn from different distributions using the parametric t-test and the emperical permutation-based method used for TACIT that uses only permutations whose test statistics have the same sign as the test statistic for the real data.
