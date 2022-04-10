Dependencies.txt lists dependencies for files in this directory beyond those listed in the main directory.

ocr_phylolm_figures.r is a modified version of ocr_phylolm.r that was used in the generation of figures. 
It includes the original ocr_phylolm.r code, for use in loading and generating data, as well as several blocks of
code below that were used to generate plots showing associations between OCR activity predictions and phenotypes.
The code included is a superset of that actually needed to generate the relavent figure panels in Kaplow, et al.
This script is not intended to be run as a whole; rather, different blocks of code were run to produce different figures.
It is provided to illustrate how figures were produced and is not intended for reuse. It also contains several hardcoded paths.

p_dist_collect.py is a Python 3 script for comparing p-value distributions for p-values associated with OCRs that are
near or not near certain genes. It takes four positional command-line arguments:
1) A list of genes of interest, one per line
2) A space-separated file with each gene near each OCR, one pair per line, with the OCRs in column 4 and the genes in column 8
3) A space-separated file with the p-value associated with each OCR (one per line), with OCRs in column 1 and p-values in column 4.
4) A file to write output, which is a CSV where the first column contains the p-values of each OCR, and the second column contains
   T if the corresponding OCR is near a gene of interest and F otherwise.
The program writes output as descibed, and additionally reports the result of a Wilcoxon rank-sum test between the set of p-values
associated with OCRs near genes of interest and the set of p-values near other genes but not genes of interests. 
OCRs near no genes (i.e. that do not occur in column 4 of input file 2) are ignored.

p_dist_collect_gene.py is a Python 3 script for comparing p-value distributions for p-values associated with genes that are
near or not near OCRS. It takes three positional command-line arguments:
1) A space-separated file with each gene near each OCR, one pair per line, with the OCRs in column 4 and the genes in column 8
2) A space-separated file with the p-value associated with  each gene (one per line), with genes in column 1 and p-values in column 4.
3) A file to write output, which is a CSV where the first column contains the p-values of each gene, and the second column contains
   T if the corresponding gene is near an OCR and F otherwise.
The program writes output as descibed, and additionally reports the result of a Wilcoxon rank-sum test between the set of p-values
associated with genes near OCRs of interest and the set of p-values associated with all other genes that occur in input file 3. 
OCRs of interest are those that appear in column 4 of input file 1.

compute_perm_pvals_kulinskaya.py is an alternate version of compute_perm_pvals.py (main directory) that, instead of using a pure one-sided test, used the method of Kulinskaya, 2008 (https://arxiv.org/abs/0810.2124), in which permulations in the opposite direction are rejected. It otherwise behaves the same as compute_perm_pvals.py, except that this script allows for multiple input directories (preceded by a count). 