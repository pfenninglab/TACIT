# Filter mouse orthologs for unique peaks
python /home/ikaplow/RegulatoryElementEvolutionProject/src/OCROrthologPrediction/src/makeFilterPeakNameScript.py --unfilteredPeakFileNameListFileName /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/Cortex_All_ATAC_out/peak/macs2/idr/optimal_set/LiftedPeaks/liftedPeakAllFileNamesConverted.txt --peakListFileName /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/Cortex_All_ATAC_out/peak/macs2/idr/optimal_set/Pfenning_bulk_Ctx_nonCDS_enhancerShort_unique.bed --peakNameCol 3 --outputFileNameSuffix unique.bed --scriptFileName /home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode/filterMouseMotorCortexOrthologsForUniqueScript.sh --codePath /home/ikaplow/RegulatoryElementEvolutionProject/src/OCROrthologPrediction/src --numFileNameElementsToRemove 2
sh /home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode/filterMouseMotorCortexOrthologsForUniqueScript.sh

# Filter macaque orthologs for peaks not overlapping mouse orthologs
grep _ma /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/Cortex_All_ATAC_out/peak/macs2/idr/optimal_set/MoMaRaBaEnhancersNRSummit.txt | cut -f1 | sed 's/_ma//g' > /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MacaqueAtac/OfM/peak/idr_reproducibility/idr.optimal_peak.inM1Loose_nonCDS_enhancerShort_nonFiltMouseOrth_peakNames.txt
python /home/ikaplow/RegulatoryElementEvolutionProject/src/OCROrthologPrediction/src/makeFilterPeakNameScript.py --unfilteredPeakFileNameListFileName /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MacaqueAtac/OfM/peak/idr_reproducibility/LiftedPeaksFixed/liftedPeakAllFileNamesConverted.txt --peakListFileName /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MacaqueAtac/OfM/peak/idr_reproducibility/idr.optimal_peak.inM1Loose_nonCDS_enhancerShort_nonFiltMouseOrth_peakNames.txt --outputFileNameSuffix nonFiltMouseOrth.bed --scriptFileName /home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode/filterMacaqueMotorCortexOrthologsForNonFiltMouseOrthScript.sh --codePath /home/ikaplow/RegulatoryElementEvolutionProject/src/OCROrthologPrediction/src --numFileNameElementsToRemove 2
sh /home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode/filterMacaqueMotorCortexOrthologsForNonFiltMouseOrthScript.sh

# Filter rat orthologs for peaks not overlapping mouse and macaque orthologs
grep _ra /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/Cortex_All_ATAC_out/peak/macs2/idr/optimal_set/MoMaRaBaEnhancersNRSummit.txt | cut -f1 | sed 's/_ra//g' > /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/RatAtac/M1_AllReps/peak/idr_reproducibility/idr.optimal_peak_nonCDS_enhancerShort_nonFiltMouseMacaqueOrth_peakNames.txt
python /home/ikaplow/RegulatoryElementEvolutionProject/src/OCROrthologPrediction/src/makeFilterPeakNameScript.py --unfilteredPeakFileNameListFileName /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/RatAtac/M1_AllReps/peak/idr_reproducibility/LiftedPeaksFixed/liftedPeakAllFileNamesConverted.txt --peakListFileName /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/RatAtac/M1_AllReps/peak/idr_reproducibility/idr.optimal_peak_nonCDS_enhancerShort_nonFiltMouseMacaqueOrth_peakNames.txt --outputFileNameSuffix nonFiltMouseMacaqueOrth.bed --scriptFileName /home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode/filterRatMotorCortexOrthologsForNonFiltMouseMacaqueOrthScript.sh --codePath /home/ikaplow/RegulatoryElementEvolutionProject/src/OCROrthologPrediction/src --numFileNameElementsToRemove 2
sh /home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode/filterRatMotorCortexOrthologsForNonFiltMouseMacaqueOrthScript.sh

# Filter bat orthologs for peaks not overlapping mouse, macaque, and rat orthologs
grep _ba /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/Cortex_All_ATAC_out/peak/macs2/idr/optimal_set/MoMaRaBaEnhancersNRSummit.txt | cut -f1 | sed 's/_ba//g' > /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/BatAtac/OfMBat1K/call-reproducibility_idr/execution/idr.optimal_peak.inM1Loose_nonCDS_enhancerShort_nonFiltMouseMacaqueRatOrth_peakNames.txt
python /home/ikaplow/RegulatoryElementEvolutionProject/src/OCROrthologPrediction/src/makeFilterPeakNameScript.py --unfilteredPeakFileNameListFileName /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/BatAtac/OfMBat1K/call-reproducibility_idr/execution/LiftedPeaksFixed/liftedPeakAllFileNamesConverted.txt --peakListFileName /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/BatAtac/OfMBat1K/call-reproducibility_idr/execution/idr.optimal_peak.inM1Loose_nonCDS_enhancerShort_nonFiltMouseMacaqueRatOrth_peakNames.txt --outputFileNameSuffix nonFiltMouseMacaqueRatOrth.bed --scriptFileName /home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode/filterBatMotorCortexOrthologsForNonFiltMouseMacaqueRatOrthScript.sh --codePath /home/ikaplow/RegulatoryElementEvolutionProject/src/OCROrthologPrediction/src --numFileNameElementsToRemove 2
sh /home/ikaplow/RegulatoryElementEvolutionProject/src/RegElEvoCode/filterBatMotorCortexOrthologsForNonFiltMouseMacaqueRatOrthScript.sh
