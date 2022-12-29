sequenceOperations.py: utilities for re-formatting sequences for convolutional neural network

MLOperations.py: utilities for evaluating machine learning models

predictNewSequences.py: evaluate predictions for a convolutional neural network on a dataset and set of labels

evaluateSingleSpeciesMotorCortexModelTestSet.sh: evaluate mouse-only motor cortex model on test set

reverseComplementFasta.py: make a fasta file that is the reverse complement of a fasta file

evaluateSequencePredictions.py: evaluate predictions for positives and negatives from a machine learning model

evaluate5SpeciesMultiSpeciesModelsTestSet.sh: evaluate multi-species motor cortex and liver models on test set

plotModelPerformanceBarGraphs5Species.m: plot the AUCs, AUNPV-Specs., and AUPRCs for a machine learning models for different evaluation criteria

makeConvertChromNamesScript.py: make a script that will convert chromosome names between naming conventions

convertChromNames.py: convert chromosome names between naming conventions

gatherPeakPredictionsAcrossSpecies.py: make an enhancer by species matrix, where each entry has the prediction for the enhancer ortholog in the species

mapCortexEnhancersAcrossZoonomia.sh: map motor cortex open chromatin regions across Zoonomia species

mapLiverEnhancersAcrossZoonomia.sh: map liver open chromatin regions across Zoonomia species

plotPredictionsVsEvolutionaryDist5Species.m: make phylogeny-matching correlations plots

species_tree.R: construct species tree using predictions at enhancer orthologs as features

ocr_heatmap_figure.r: plot heatmap of OCR activity predictions

filterMotorCortexPeakOrthologs.sh: filter OCR orthologs so that no regions are repeated

makeSummitCenteredFastaFromNarrowPeak.py: convert a narrowPeak file to a summit-centered fasta file

makeSummitCenteredFastaFromNarrowPeakScript.py: make a script that runs makeSummitCenteredFastaFromNarrowPeak.py on a list of narrowPeak files

makeFIMOListScript.py: make a script that runs FIMO from the MEME suite on a list of fasta files

logistic_train_modified.py: train, evaluate, or make predictions with a logistic regression model with motif score features

makeLogisticTrainModifiedScript.py: make a script that uses logistic_train_modified.py to make logistic model predictions for features from a list of FIMO files
