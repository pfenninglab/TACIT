# TACIT Enhancer-Phenotype Association
Enhancer-Phenotype Association from Tissue-Aware Conservation Inference Toolkit.

We used this with open chromatin regions (OCRs) / predictions of open chromatin status, but the pipeline is also amenable to predictions of other proxies for enhancer activity. 

This pipeline has been tested on Centos versions 7 and 8. However, we expect it to work on other distributions of Linux, MacOS, and Windows.

## Programs for OCR-Phenotype association with phylo(g)lm:
`ocr_phylolm.r`:    Performs OCR-phenotype assosiations for continuous traits, using phylolm, possibly with permulations. See the header comment of this program for usage details.
                  
`ocr_phyloglm.r`:   Performs OCR-phenotype assosiations for binary traits, using phyloglm, possibly with permulations. See the header comment of this program for usage details. 

`fast_bin_perm.r`: Our fast implementation of binary permulations (1). It contains the function `fastSimBinPhenoVec`, which is called by `ocr_phyloglm.r` and can also be used to generate permulations of binary traits for other statistical tests. The function `makeLeafMap` precomputes the required data structure. The other functions are called only internally.


## Utilities:
`compute_perm_pvals.py`: Script for computing permulation p-values using the output from `ocr_phylo(g)lm.r`. See the header comment of this script for usage details.

`compute_perm_pvals_kulinskaya.py` is an alternate version of compute_perm_pvals.py that, instead of using a pure one-sided test, used the method of Kulinskaya, 2008 (https://arxiv.org/abs/0810.2124), in which permulations in the opposite direction are rejected. It otherwise behaves the same as compute_perm_pvals.py, except that this script allows for multiple input directories (preceded by a count).  We recommend using this method instead of compute_perm_pvals.py for traits that can be easily permulated, as it is more principled but slower due to its rejection of permulations in the opposite direction.


## Dependencies:
### For R programs:
R >= 4.0

ape package >= 5.4.1 (2)

phylolm package >= 2.6.2 (3)

geiger package >= 2.0.7 (4)

Rcpp package >= 1.0.6 (5)

stringr package >= 1.4.0

### For Python programs:
Python 3

## Scripts for generating figures in Kaplow, Lawler, Schaffer, _et al_.:

evaluationScriptsRERconverge: scripts for preparing data for and running RERconverge

evaluationScriptsCortexLiverModels: scripts for evaluating and making predictions with motor cortex and liver models

evaluationScriptsPVModels: scripts for evaluating and making predictions with PV models

evaluationScriptsRetinaModels: scripts for evaluating and making predictions with retina models

evaluationScriptsEnhancerPhenotypeAssociations: scripts for evaluating and analyzing open chromatin region-phenotype associations

evaluationScriptsHiC: scripts for evaluating long-range interactions between phenotype-associated open chromatin regions and phenotype-relevant genes

Dependencies specific to those scripts are listed in the individual directories.  Models can be found at http://daphne.compbio.cs.cmu.edu/files/ikaplow/TACITSupplement/.

## Contact
Daniel Schaffer (dschaffe@andrew.cmu.edu)

Irene Kaplow (ikaplow@cs.cmu.edu)

Alyssa Lawler (alawler@andrew.cmu.edu)

Chaitanya Srinivasan (csriniv1@andrew.cmu.edu)

Carson Sestili (csestili@andrew.cmu.edu)

Andreas Pfenning (apfenning@cmu.edu)

## References
1. I. M. Kaplow, A. J. Lawler, D. E. Schaffer, C. Srinivasan, M. E. Wirthlin, B. N. Phan, X. Zhang, K. Foley, A. R. Brown, Zoonomia Consortium, W. K. Meyer, A. R. Pfenning, Relating enhancer genetic variation across mammals to complex phenotypes using machine learning. *bioRxiv*. https://www.biorxiv.org/content/10.1101/2022.08.26.505436v1 (2022).
2. D. Eddelbuettel, R. François, J. Allaire, K. Ushey, Q. Kou, N. Russel, J. Chambers, D. Bates, Rcpp: Seamless R and C++ integration. *J. Stat. Softw*. **40**, 1–18 (2011).
3. L. S. T. Ho, C. Ané, A linear-time algorithm for Gaussian and non-Gaussian trait evolution models. *Syst. Biol*. **63**, 397–408 (2014).
4. E. Kulinskaya, On two-sided p-values for non-symmetric distributions. *arXiv*. https://arxiv.org/abs/0810.2124 (2008).
5. E. Paradis, K. Schliep, ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. *Bioinformatics*. **35**, 526–528 (2019).
6. M. W. Pennell, J. M. Eastman, G. J. Slater, J. W. Brown, J. C. Uyeda, R. G. FitzJohn, M. E. Alfaro, L. J. Harmon, geiger v2.0: an expanded suite of methods for fitting macroevolutionary models to phylogenetic trees. *Bioinformatics*. **30** (2014), pp. 2216–2218.
7. E. Saputra, A. Kowalczyk, L. Cusick, N. Clark, M. Chikina, Phylogenetic Permulations: A Statistically Rigorous Approach to Measure Confidence in Associations in a Phylogenetic Context. *Mol. Biol. Evol*. **38**, 3004–3021 (2021).
