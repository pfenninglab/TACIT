# TACIT OCR-Phenotype Association
OCR-Phenotype Association from Tissue-Aware Conservation Inference Toolkit

## Programs for OCR-Phenotype association with phylo(g)lm:
ocr_phylolm.r:    Performs OCR-phenotype assosiations for continuous traits, using phylolm, possibly with permulations.
                  See the header comment of this program for usage details.
ocr_phyloglm.r:   Performs OCR-phenotype assosiations for binary traits, using phyloglm, possibly with permulations.
                  See the header comment of this program for usage details.     
fast_bin_perm.r: Our fast implementation of binary permulations, for use with ocr_phyloglm.r. See the method signatures
                 in thie file for usage details; it can also be used with other statistical tests. 

## Utilities:
compute_perm_pvals.py: Script for computing permulation p-values using the output from ocr_phylo(g)lm.r.
                       See the header comment of this script for usage details.

## Dependencies:
### For R programs:
R >= 4.0
ape package >= 5.4.1
phylolm package >= 2.6.2
geiger package >= 2.0.7
Rcpp package >= 1.0.6
stringr package >= 1.4.0

### For Python programs:
Python 3

