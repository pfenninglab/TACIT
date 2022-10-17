import os
import pickle

import matplotlib.pyplot as plt
import numpy as np
import scipy
from tqdm import tqdm

import visualization
import models


# Global constants

REDUCER_TYPE = 'pca' # 'pca' or 'umap'
if REDUCER_TYPE not in ['pca', 'umap']:
	raise NotImplementedError()
PLOT_EXTENSION = 'both' # 'both', 'png', or 'svg'

FIG_WIDTH_CM = 4.6231


# Model-specific constants

# PV models

# # multispecies PV model
# # MODEL_PATH = '/projects/pfenninggroup/mouseCxStr/NeuronSubtypeATAC/Zoonomia_CNN/multispecies_PV/models/FINAL_modelmultiPVi.h5'
# # LAYER_NAME = 'activation_1'
# FIT_DATA_KEYS = ['combined_pos_val', 'combined_neg_val']
# ACTIVATIONS_DIR = '/home/csestili/data/tacit_viz/pv_multispecies/'

# # # mouse-only PV model
# # MODEL_PATH = '/projects/pfenninggroup/mouseCxStr/NeuronSubtypeATAC/Zoonomia_CNN/mouse_PV/models/FINAL_modelPV3e.h5'
# # # [...] -> maxpool -> flatten -> dense(300) -> activation_9 -> [...]
# # LAYER_NAME = 'activation_9'
# # FIT_DATA_KEYS = ['mouse_pos_val', 'mouse_neg_val']
# # ACTIVATIONS_DIR = '/home/csestili/data/tacit_viz/pv_mouseonly/'

# CAN_LOAD_MODEL = True

# # PV_DATA_DIR = '/projects/pfenninggroup/mouseCxStr/NeuronSubtypeATAC/Zoonomia_CNN/multispecies_PV/FinalModelData'
# # PV_EVAL_DIR = '/projects/pfenninggroup/mouseCxStr/NeuronSubtypeATAC/Zoonomia_CNN/evaluations/PV/'
# # DATA_MAPPING = {
# # 	'mouse_neg_val': os.path.join(PV_DATA_DIR, 'mouse_PV_neg_VAL.fa'),
# # 	'combined_pos_val': os.path.join(PV_DATA_DIR, 'combined_PV_pos_VAL.fa'),
# # 	'combined_neg_val': os.path.join(PV_DATA_DIR, 'combined_PV_neg_VAL.fa'),
# # 	'human_pos_train': os.path.join(PV_DATA_DIR, 'human_PV_pos_TRAIN.fa'),
# # 	'human_pos_val': os.path.join(PV_DATA_DIR, 'human_PV_pos_VAL.fa'),
# # 	'human_pos_test': os.path.join(PV_DATA_DIR, 'human_PV_pos_TEST.fa'),
# # 	'mouse_pos_train': os.path.join(PV_DATA_DIR, 'mouse_PV_pos_TRAIN.fa'),
# # 	'mouse_pos_val': os.path.join(PV_DATA_DIR, 'mouse_PV_pos_VAL.fa'),
# # 	'mouse_pos_test': os.path.join(PV_DATA_DIR, 'mouse_PV_pos_TEST.fa'),
# # 	'human_neg_neoe_all': os.path.join(PV_EVAL_DIR, 'Eval4_mm10.fa'),
# # 	'mouse_neg_neoe_all': os.path.join(PV_EVAL_DIR, 'Eval2_hg38.fa'),
# # 	'dolphin': '/projects/pfenninggroup/mouseCxStr/NeuronSubtypeATAC/Zoonomia_CNN/predictions/mousePV_per_species/fasta/mouseReproduciblePV_mm10_Tursiops_truncatus.fa',
# # 	'rabbit': '/projects/pfenninggroup/mouseCxStr/NeuronSubtypeATAC/Zoonomia_CNN/predictions/mousePV_per_species/fasta/glires/mouseReproduciblePV_mm10_Oryctolagus_cuniculus.fa'
# # }

# VISUALIZATION_MAPPING = {
# 	'fit_data': FIT_DATA_KEYS,
# 	'plot_data': [
# 		# {
# 		# 	'title': 'Positives',
# 		# 	'groups': [
# 		# 		{'id': 4, 'name': 'Human +', 'sets': ['human_pos_train', 'human_pos_val', 'human_pos_test']},
# 		# 		{'id': 5, 'name': 'Mouse +', 'sets': ['mouse_pos_train', 'mouse_pos_val', 'mouse_pos_test']}
# 		# 	],
# 		# 	'add_histogram': True
# 		# },
# 		# {
# 		# 	'title': 'Human',
# 		# 	'groups': [
# 		# 		{'id': 0, 'name': 'Human - (neoe)', 'sets': ['human_neg_neoe_all']},
# 		# 		{'id': 4, 'name': 'Human +', 'sets': ['human_pos_train', 'human_pos_val', 'human_pos_test']}
# 		# 	],
# 		# 	'add_histogram': True
# 		# },
# 		# {
# 		# 	'title': 'Mouse',
# 		# 	'groups': [
# 		# 		{'id': 1, 'name': 'Mouse - (neoe)', 'sets': ['mouse_neg_neoe_all']},
# 		# 		{'id': 5, 'name': 'Mouse +', 'sets': ['mouse_pos_train', 'mouse_pos_val', 'mouse_pos_test']}
# 		# 	],
# 		# 	'add_histogram': True
# 		# },
# 		# # plot order 1
# 		# {
# 		# 	'title': 'Mammals',
# 		# 	'groups': [
# 		# 		{'id': 5, 'name': 'Mouse +', 'sets': ['mouse_pos_train', 'mouse_pos_val', 'mouse_pos_test']},
# 		# 		{'id': 4, 'name': 'Human +', 'sets': ['human_pos_train', 'human_pos_val', 'human_pos_test']},
# 		# 		{'id': 3, 'name': 'Rabbit', 'sets': ['rabbit']},
# 		# 		{'id': 2, 'name': 'Dolphin', 'sets': ['dolphin']},
# 		# 		{'id': 1, 'name': 'Mouse - (neoe)', 'sets': ['mouse_neg_neoe_all']},
# 		# 		{'id': 0, 'name': 'Human - (neoe)', 'sets': ['human_neg_neoe_all']}
# 		# 	],
# 		# 	'add_violinplot': True,
# 		# 	'add_ranksum_table': True
# 		# },
# 		# plot order 2
# 		{
# 			'title': 'Mammals (order 2)',
# 			'groups': [
# 				{'id': 5, 'name': 'Mouse +', 'sets': ['mouse_pos_train', 'mouse_pos_val', 'mouse_pos_test']},
# 				{'id': 4, 'name': 'Mouse - (neoe)', 'sets': ['mouse_neg_neoe_all']},
# 				{'id': 3, 'name': 'Human +', 'sets': ['human_pos_train', 'human_pos_val', 'human_pos_test']},
# 				{'id': 2, 'name': 'Human - (neoe)', 'sets': ['human_neg_neoe_all']},
# 				{'id': 1, 'name': 'Rabbit', 'sets': ['rabbit']},
# 				{'id': 0, 'name': 'Dolphin', 'sets': ['dolphin']}
# 			],
# 			'add_violinplot': True,
# 			'add_ranksum_table': True
# 		}
# 	]
# }

# Motor Cortex models

# CAN_LOAD_MODEL = False
# # placeholder model
# MODEL_PATH = '/projects/pfenninggroup/mouseCxStr/NeuronSubtypeATAC/Zoonomia_CNN/multispecies_PV/models/FINAL_modelmultiPVi.h5'
# LAYER_NAME = 'activation_1'

# # multispecies MC model
# #MODEL_PATH = '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/cortexEnhancerVsOtherSpeciesEnhLoosePeakLessThan1kb500bp_data_train.KerasModels/hybrid_cortexEnhancerShort_enhLooseNeg_500bp_conv5MuchMuchMoreFiltNoInterPoolTwoDenseLargeDenseVeryHighMomLowDropL2SmallBatchPretrainBal.hdf5'
# FIT_DATA_KEYS = ['mouse_pos_val', 'mouse_neg_val',
# 	'mac_pos_val', 'rat_pos_val', 'bat_pos_val',
# 	'mac_neg_val', 'rat_neg_val', 'bat_neg_val']
# ACTIVATIONS_DIR = '/home/csestili/data/tacit_viz/mc_multispecies/'

# # # mouse-only MC model
# # #MODEL_PATH = '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/cortexEnhancerVsOtherSpeciesEnhMouseMacaqueRatBatPeakLessThan1kb_data_train.KerasModels/cortexEnhancerShort_multiSpecies_enhLooseNeg_500bp_conv5TotalTotalFiltNoInterPoolVeryHighMomL2SmallBatchPretrainBal.hdf5'

# # FIT_DATA_KEYS = ['mouse_pos_val', 'mouse_neg_val']
# # ACTIVATIONS_DIR = '/home/csestili/data/tacit_viz/mc_mouseonly/'

# DATA_MAPPING = {
# 	'mouse_pos_val': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/Cortex_All_ATAC_out/peak/macs2/idr/optimal_set/Pfenning_bulk_Ctx_nonCDS_enhancerShort_valid_summitPlusMinus250bp.fa',
# 	'mouse_neg_val': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MacaqueAtac/OfM/peak/idr_reproducibility/idr.optimal_peak.inM1Loose_nonCDS_enhancerShort_GenBankNames_mm10_summitExtendedMin50Max2XProtect5_nonMouseCortex_valid_andRat_andBat_summitPlusMinus250bp.fa',
# 	'mouse_pos': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/Cortex_All_ATAC_out/peak/macs2/idr/optimal_set/Pfenning_bulk_Ctx_nonCDS_enhancerShort_unique_summitPlusMinus250bp.fa',
# 	'mac_pos_train': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MacaqueAtac/OfM/peak/idr_reproducibility/idr.optimal_peak.inM1_nonCDS_enhancerShort_train_summitPlusMinus250bp.fa', 
# 	'mac_pos_val': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MacaqueAtac/OfM/peak/idr_reproducibility/idr.optimal_peak.inM1_nonCDS_enhancerShort_valid_summitPlusMinus250bp.fa',
# 	'mac_pos_test': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MacaqueAtac/OfM/peak/idr_reproducibility/idr.optimal_peak.inM1_nonCDS_enhancerShort_test_summitPlusMinus250bp.fa',
# 	'rat_pos_train': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/RatAtac/M1_AllReps/peak/idr_reproducibility/idr.optimal_peak_nonCDS_enhancerShort_train_summitPlusMinus250bp.fa',
# 	'rat_pos_val': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/RatAtac/M1_AllReps/peak/idr_reproducibility/idr.optimal_peak_nonCDS_enhancerShort_valid_summitPlusMinus250bp.fa',
# 	'rat_pos_test': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/RatAtac/M1_AllReps/peak/idr_reproducibility/idr.optimal_peak_nonCDS_enhancerShort_test_summitPlusMinus250bp.fa',
# 	'bat_pos_train': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/BatAtac/OfMBat1K/call-reproducibility_idr/execution/idr.optimal_peak.inM1_nonCDS_enhancerShort_train_summitPlusMinus250bp.fa',
# 	'bat_pos_val': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/BatAtac/OfMBat1K/call-reproducibility_idr/execution/idr.optimal_peak.inM1_nonCDS_enhancerShort_valid_summitPlusMinus250bp.fa',
# 	'bat_pos_test': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/BatAtac/OfMBat1K/call-reproducibility_idr/execution/idr.optimal_peak.inM1_nonCDS_enhancerShort_test_summitPlusMinus250bp.fa',
# 	'mouse_neg': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MacaqueAtac/OfM/peak/idr_reproducibility/idr.optimal_peak.inM1Loose_nonCDS_enhancerShort_GenBankNames_mm10_summitExtendedMin50Max2XProtect5_nonMouseCortex_andRat_andBat.plusMinus250bp.fa',
# 	'mac_neg_train': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/Cortex_All_ATAC_out/peak/macs2/idr/optimal_set/Pfenning_bulk_Ctx_nonCDS_enhancerShort_rheMac8_summitExtendedMin50Max2XProtect5_UCSCNames_nonMacaqueCortex_train_andRat_andBat_summitPlusMinus250bp.fa',
# 	'mac_neg_val': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/Cortex_All_ATAC_out/peak/macs2/idr/optimal_set/Pfenning_bulk_Ctx_nonCDS_enhancerShort_rheMac8_summitExtendedMin50Max2XProtect5_UCSCNames_nonMacaqueCortex_valid_andRat_andBat_summitPlusMinus250bp.fa',
# 	'mac_neg_test': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/Cortex_All_ATAC_out/peak/macs2/idr/optimal_set/Pfenning_bulk_Ctx_nonCDS_enhancerShort_rheMac8_summitExtendedMin50Max2XProtect5_UCSCNames_nonMacaqueCortex_test_andRat_andBat_summitPlusMinus250bp.fa',
# 	'rat_neg_train': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/Cortex_All_ATAC_out/peak/macs2/idr/optimal_set/Pfenning_bulk_Ctx_nonCDS_enhancerShort_rn6_summitExtendedMin50Max2XProtect5_nonRatCortex_train_andMacaque_andBat_summitPlusMinus250bp.fa',
# 	'rat_neg_val': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/Cortex_All_ATAC_out/peak/macs2/idr/optimal_set/Pfenning_bulk_Ctx_nonCDS_enhancerShort_rn6_summitExtendedMin50Max2XProtect5_nonRatCortex_valid_andMacaque_andBat_summitPlusMinus250bp.fa',
# 	'rat_neg_test': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/Cortex_All_ATAC_out/peak/macs2/idr/optimal_set/Pfenning_bulk_Ctx_nonCDS_enhancerShort_rn6_summitExtendedMin50Max2XProtect5_nonRatCortex_test_andMacaque_andBat_summitPlusMinus250bp.fa',
# 	'bat_neg_train': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/Cortex_All_ATAC_out/peak/macs2/idr/optimal_set/Pfenning_bulk_Ctx_nonCDS_enhancerShort_200MBat_summitExtendedMin50Max2XProtect5_HLrouAeg4_RefSeqNames_nonBatCortex_train_andMacaque_andRat_summitPlusMinus250bp.fa',
# 	'bat_neg_val': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/Cortex_All_ATAC_out/peak/macs2/idr/optimal_set/Pfenning_bulk_Ctx_nonCDS_enhancerShort_200MBat_summitExtendedMin50Max2XProtect5_HLrouAeg4_RefSeqNames_nonBatCortex_valid_andMacaque_andRat_summitPlusMinus250bp.fa',
# 	'bat_neg_test': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/Cortex_All_ATAC_out/peak/macs2/idr/optimal_set/Pfenning_bulk_Ctx_nonCDS_enhancerShort_200MBat_summitExtendedMin50Max2XProtect5_HLrouAeg4_RefSeqNames_nonBatCortex_test_andMacaque_andRat_summitPlusMinus250bp.fa',
# 	'rabbit': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/Cortex_All_ATAC_out/peak/macs2/idr/optimal_set/LiftedPeaks/Pfenning_bulk_Ctx_nonCDS_enhancerShort_Oryctolagus_cuniculus_summitExtendedMin50Max2XProtect5_RefSeqNames.plusMinus250bp.fa',
# 	'dolphin': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/Cortex_All_ATAC_out/peak/macs2/idr/optimal_set/LiftedPeaks/Pfenning_bulk_Ctx_nonCDS_enhancerShort_Tursiops_truncatus_summitExtendedMin50Max2XProtect5_RefSeqNames.plusMinus250bp.fa'
# }

# VISUALIZATION_MAPPING = {
# 	'fit_data': FIT_DATA_KEYS,
# 	'plot_data': [
# 		# {
# 		# 	'title': 'Positives',
# 		# 	'groups': [
# 		# 		{'id': 4, 'name': 'Human +', 'sets': ['human_pos_train', 'human_pos_val', 'human_pos_test']},
# 		# 		{'id': 5, 'name': 'Mouse +', 'sets': ['mouse_pos_train', 'mouse_pos_val', 'mouse_pos_test']}
# 		# 	],
# 		# 	'add_histogram': True
# 		# },
# 		# {
# 		# 	'title': 'Human',
# 		# 	'groups': [
# 		# 		{'id': 0, 'name': 'Human - (neoe)', 'sets': ['human_neg_neoe_all']},
# 		# 		{'id': 4, 'name': 'Human +', 'sets': ['human_pos_train', 'human_pos_val', 'human_pos_test']}
# 		# 	],
# 		# 	'add_histogram': True
# 		# },
# 		# {
# 		# 	'title': 'Mouse',
# 		# 	'groups': [
# 		# 		{'id': 1, 'name': 'Mouse - (neoe)', 'sets': ['mouse_neg_neoe_all']},
# 		# 		{'id': 5, 'name': 'Mouse +', 'sets': ['mouse_pos_train', 'mouse_pos_val', 'mouse_pos_test']}
# 		# 	],
# 		# 	'add_histogram': True
# 		# },
# 		{
# 			'title': 'Mammals',
# 			'groups': [
# 				{'id': 9, 'name': 'Mouse +', 'sets': ['mouse_pos']},
# 				{'id': 8, 'name': 'Mouse -', 'sets': ['mouse_neg']},
# 				{'id': 7, 'name': 'Rat +', 'sets': ['rat_pos_train', 'rat_pos_val', 'rat_pos_test']},
# 				{'id': 6, 'name': 'Rat -', 'sets': ['rat_neg_train', 'rat_neg_val', 'rat_neg_test']},
# 				{'id': 5, 'name': 'Macaque +', 'sets': ['mac_pos_train', 'mac_pos_val', 'mac_pos_test']},
# 				{'id': 4, 'name': 'Macaque -', 'sets': ['mac_neg_train', 'mac_neg_val', 'mac_neg_test']},
# 				{'id': 3, 'name': 'Bat +', 'sets': ['bat_pos_train', 'bat_pos_val', 'bat_pos_test']},
# 				{'id': 2, 'name': 'Bat -', 'sets': ['bat_neg_train', 'bat_neg_val', 'bat_neg_test']},
# 				{'id': 1, 'name': 'Rabbit', 'sets': ['rabbit']},
# 				{'id': 0, 'name': 'Dolphin', 'sets': ['dolphin']}
# 			],
# 			'add_violinplot': True,
# 			'add_ranksum_table': True
# 		}
# 	]
# }



# liver model

# # real model path
# # MODEL_PATH = '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/liverEnhancerVsOtherSpeciesEnhMouseMacaqueRatCowPigPeakLessThan1kb_data_train.KerasModels/liver_mouseMacaqueRatCowPig_otherSpeciesEnhLooseNeg_500bp_conv5AllFiltNoInterPoolVeryHighMomL2SmallBatchPretrainBal_model.hdf5'

# # placeholder model
# MODEL_PATH = '/projects/pfenninggroup/mouseCxStr/NeuronSubtypeATAC/Zoonomia_CNN/multispecies_PV/models/FINAL_modelmultiPVi.h5'
# LAYER_NAME = 'activation_1'
# CAN_LOAD_MODEL = False

# FIT_DATA_KEYS = ['mouse_pos_val', 'mouse_neg_val',
# 	'mac_pos_val', 'rat_pos_val', 'cow_pos_val', 'pig_pos_val',
# 	'mac_neg_val', 'rat_neg_val', 'cow_neg_val', 'pig_neg_val']
# ACTIVATIONS_DIR = '/home/csestili/data/tacit_viz/liver_multispecies/'

# DATA_MAPPING = {
# 	'mouse_pos_val': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/LiverPfenning_ATAC_out/atac/68eb27d1-23c7-4d42-80f9-b3dd46f5e18c/call-reproducibility_idr/execution/idr.optimal_peak.inLiuAll_nonCDS_enhancerShort_valid_summitPlusMinus250bp.fa',
# 	'mouse_neg_val': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MacaqueAtac/Liver/cromwell-executions/atac/1f5fd9f0-d389-4184-99cb-c1125fd7f064/call-reproducibility_idr/execution/optimal_peak_nonCDS_enhancerShort_mm10Fixed_summitExtendedMin50Max2Protect5_nonMouseLiver_andRat_andCow_andPig_valid_summitPlusMinus250bp.fa',
# 	'mouse_pos': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/LiverPfenning_ATAC_out/atac/68eb27d1-23c7-4d42-80f9-b3dd46f5e18c/call-reproducibility_idr/execution/idr.optimal_peak.narrowPeak_inLiuAll_nonCDS_enhancerShort_summitPlusMinus250bp.fa',
# 	'mac_pos_train': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MacaqueAtac/Liver/cromwell-executions/atac/1f5fd9f0-d389-4184-99cb-c1125fd7f064/call-reproducibility_idr/execution/optimal_peak_nonCDS_enhancerShort_train_summitPlusMinus250bp.fa', 
# 	'mac_pos_val': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MacaqueAtac/Liver/cromwell-executions/atac/1f5fd9f0-d389-4184-99cb-c1125fd7f064/call-reproducibility_idr/execution/optimal_peak_nonCDS_enhancerShort_valid_summitPlusMinus250bp.fa',
# 	'mac_pos_test': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MacaqueAtac/Liver/cromwell-executions/atac/1f5fd9f0-d389-4184-99cb-c1125fd7f064/call-reproducibility_idr/execution/optimal_peak_nonCDS_enhancerShort_test_summitPlusMinus250bp.fa',
# 	'rat_pos_train': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/RatAtac/Liver_GoodReps/peak/idr_reproducibility/idr.optimal_peak_nonCDS_enhancerShort_train_summitPlusMinus250bp.fa',
# 	'rat_pos_val': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/RatAtac/Liver_GoodReps/peak/idr_reproducibility/idr.optimal_peak_nonCDS_enhancerShort_valid_summitPlusMinus250bp.fa',
# 	'rat_pos_test': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/RatAtac/Liver_GoodReps/peak/idr_reproducibility/idr.optimal_peak_nonCDS_enhancerShort_test_summitPlusMinus250bp.fa',
# 	'cow_pos_train': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/CowDNase/Liver/atac/888ae501-7469-4afa-8890-e36f14f6c83b/peak/idr_reproducibility/idr.optimal_peak_nonCDS_enhancerShort_train_summitPlusMinus250bp.fa',
# 	'cow_pos_val': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/CowDNase/Liver/atac/888ae501-7469-4afa-8890-e36f14f6c83b/peak/idr_reproducibility/idr.optimal_peak_nonCDS_enhancerShort_valid_summitPlusMinus250bp.fa',
# 	'cow_pos_test': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/CowDNase/Liver/atac/888ae501-7469-4afa-8890-e36f14f6c83b/peak/idr_reproducibility/idr.optimal_peak_nonCDS_enhancerShort_test_summitPlusMinus250bp.fa',
# 	'pig_pos_train': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/PigDNase/Liver/atac/568838bc-4846-4491-8787-45b51e65bfee/peak/idr_reproducibility/idr.optimal_peak_nonCDS_enhancerShort_train_summitPlusMinus250bp.fa',
# 	'pig_pos_val': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/PigDNase/Liver/atac/568838bc-4846-4491-8787-45b51e65bfee/peak/idr_reproducibility/idr.optimal_peak_nonCDS_enhancerShort_valid_summitPlusMinus250bp.fa',
# 	'pig_pos_test': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/PigDNase/Liver/atac/568838bc-4846-4491-8787-45b51e65bfee/peak/idr_reproducibility/idr.optimal_peak_nonCDS_enhancerShort_test_summitPlusMinus250bp.fa',
# 	'mouse_neg': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MacaqueAtac/Liver/cromwell-executions/atac/1f5fd9f0-d389-4184-99cb-c1125fd7f064/call-reproducibility_idr/execution/optimal_peak_nonCDS_enhancerShort_mm10Fixed_summitExtendedMin50Max2Protect5_nonMouseLiver_andRat_andCow_andPig.plusMinus250bp.fa',
# 	'mac_neg_train': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/LiverPfenning_ATAC_out/atac/68eb27d1-23c7-4d42-80f9-b3dd46f5e18c/call-reproducibility_idr/execution/idr.optimal_peak.inLiuAllLoose_nonCDS_enhancerShort_rheMac8Fixed_summitExtendedMin50Max2XProtect5_UCSCNames_nonMacaqueLiver_train_andRat_andCow_andPig_summitPlusMinus250bp.fa',
# 	'mac_neg_val': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/LiverPfenning_ATAC_out/atac/68eb27d1-23c7-4d42-80f9-b3dd46f5e18c/call-reproducibility_idr/execution/idr.optimal_peak.inLiuAllLoose_nonCDS_enhancerShort_rheMac8Fixed_summitExtendedMin50Max2XProtect5_UCSCNames_nonMacaqueLiver_valid_andRat_andCow_andPig_summitPlusMinus250bp.fa',
# 	'mac_neg_test': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/LiverPfenning_ATAC_out/atac/68eb27d1-23c7-4d42-80f9-b3dd46f5e18c/call-reproducibility_idr/execution/idr.optimal_peak.inLiuAllLoose_nonCDS_enhancerShort_rheMac8Fixed_summitExtendedMin50Max2XProtect5_UCSCNames_nonMacaqueLiver_test_andRat_andCow_andPig_summitPlusMinus250bp.fa',
# 	'rat_neg_train': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/LiverPfenning_ATAC_out/atac/68eb27d1-23c7-4d42-80f9-b3dd46f5e18c/call-reproducibility_idr/execution/idr.optimal_peak.inLiuAllLoose_nonCDS_enhancerShort_rn6Fixed_summitExtendedMin50Max2XProtect5_nonRatLiver_train_andMacaque_andCow_andPig_summitPlusMinus250bp.fa',
# 	'rat_neg_val': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/LiverPfenning_ATAC_out/atac/68eb27d1-23c7-4d42-80f9-b3dd46f5e18c/call-reproducibility_idr/execution/idr.optimal_peak.inLiuAllLoose_nonCDS_enhancerShort_rn6Fixed_summitExtendedMin50Max2XProtect5_nonRatLiver_valid_andMacaque_andCow_andPig_summitPlusMinus250bp.fa',
# 	'rat_neg_test': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/LiverPfenning_ATAC_out/atac/68eb27d1-23c7-4d42-80f9-b3dd46f5e18c/call-reproducibility_idr/execution/idr.optimal_peak.inLiuAllLoose_nonCDS_enhancerShort_rn6Fixed_summitExtendedMin50Max2XProtect5_nonRatLiver_test_andMacaque_andCow_andPig_summitPlusMinus250bp.fa',
# 	'cow_neg_train': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/LiverPfenning_ATAC_out/atac/68eb27d1-23c7-4d42-80f9-b3dd46f5e18c/call-reproducibility_idr/execution/idr.optimal_peak.inLiuAllLoose_nonCDS_enhancerShort_bosTau2015Fixed_summitExtendedMin50Max2XProtect5_RefSeqNames_nonCowLiver_train_andMacaque_andRat_andPig_summitPlusMinus250bp.fa',
# 	'cow_neg_val': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/LiverPfenning_ATAC_out/atac/68eb27d1-23c7-4d42-80f9-b3dd46f5e18c/call-reproducibility_idr/execution/idr.optimal_peak.inLiuAllLoose_nonCDS_enhancerShort_bosTau2015Fixed_summitExtendedMin50Max2XProtect5_RefSeqNames_nonCowLiver_valid_andMacaque_andRat_andPig_summitPlusMinus250bp.fa',
# 	'cow_neg_test': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/LiverPfenning_ATAC_out/atac/68eb27d1-23c7-4d42-80f9-b3dd46f5e18c/call-reproducibility_idr/execution/idr.optimal_peak.inLiuAllLoose_nonCDS_enhancerShort_bosTau2015Fixed_summitExtendedMin50Max2XProtect5_RefSeqNames_nonCowLiver_test_andMacaque_andRat_andPig_summitPlusMinus250bp.fa',
# 	'pig_neg_train': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/LiverPfenning_ATAC_out/atac/68eb27d1-23c7-4d42-80f9-b3dd46f5e18c/call-reproducibility_idr/execution/idr.optimal_peak.inLiuAllLoose_nonCDS_enhancerShort_susScr3Fixed_summitExtendedMin50Max2XProtect5_UCSCNames_nonPigLiver_train_andMacaque_andRat_andCow_summitPlusMinus250bp.fa',
# 	'pig_neg_val': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/LiverPfenning_ATAC_out/atac/68eb27d1-23c7-4d42-80f9-b3dd46f5e18c/call-reproducibility_idr/execution/idr.optimal_peak.inLiuAllLoose_nonCDS_enhancerShort_susScr3Fixed_summitExtendedMin50Max2XProtect5_UCSCNames_nonPigLiver_valid_andMacaque_andRat_andCow_summitPlusMinus250bp.fa',
# 	'pig_neg_test': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/LiverPfenning_ATAC_out/atac/68eb27d1-23c7-4d42-80f9-b3dd46f5e18c/call-reproducibility_idr/execution/idr.optimal_peak.inLiuAllLoose_nonCDS_enhancerShort_susScr3Fixed_summitExtendedMin50Max2XProtect5_UCSCNames_nonPigLiver_test_andMacaque_andRat_andCow_summitPlusMinus250bp.fa',
# 	'rabbit': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/LiverPfenning_ATAC_out/atac/68eb27d1-23c7-4d42-80f9-b3dd46f5e18c/call-reproducibility_idr/execution/LiftedPeaksFixed/idr.optimal_peak.narrowPeak_inLiuAll_nonCDS_enhancerShort_Oryctolagus_cuniculus_summitExtendedMin50Max2XProtect5_RefSeqNames.plusMinus250bp.fa',
# 	'dolphin': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/LiverPfenning_ATAC_out/atac/68eb27d1-23c7-4d42-80f9-b3dd46f5e18c/call-reproducibility_idr/execution/LiftedPeaksFixed/idr.optimal_peak.narrowPeak_inLiuAll_nonCDS_enhancerShort_Tursiops_truncatus_summitExtendedMin50Max2XProtect5_RefSeqNames.plusMinus250bp.fa'
# }

# VISUALIZATION_MAPPING = {
# 	'fit_data': FIT_DATA_KEYS,
# 	'plot_data': [
# 		{
# 			'title': 'Mammals',
# 			'groups': [
# 				{'id': 11, 'name': 'Mouse +', 'sets': ['mouse_pos']},
# 				{'id': 10, 'name': 'Mouse -', 'sets': ['mouse_neg']},
# 				{'id': 9, 'name': 'Rat +', 'sets': ['rat_pos_train', 'rat_pos_val', 'rat_pos_test']},
# 				{'id': 8, 'name': 'Rat -', 'sets': ['rat_neg_train', 'rat_neg_val', 'rat_neg_test']},
# 				{'id': 7, 'name': 'Macaque +', 'sets': ['mac_pos_train', 'mac_pos_val', 'mac_pos_test']},
# 				{'id': 6, 'name': 'Macaque -', 'sets': ['mac_neg_train', 'mac_neg_val', 'mac_neg_test']},
# 				{'id': 5, 'name': 'Cow +', 'sets': ['cow_pos_train', 'cow_pos_val', 'cow_pos_test']},
# 				{'id': 4, 'name': 'Cow -', 'sets': ['cow_neg_train', 'cow_neg_val', 'cow_neg_test']},
# 				{'id': 3, 'name': 'Pig +', 'sets': ['pig_pos_train', 'pig_pos_val', 'pig_pos_test']},
# 				{'id': 2, 'name': 'Pig -', 'sets': ['pig_neg_train', 'pig_neg_val', 'pig_neg_test']},
# 				{'id': 1, 'name': 'Rabbit', 'sets': ['rabbit']},
# 				{'id': 0, 'name': 'Dolphin', 'sets': ['dolphin']}
# 			],
# 			'add_violinplot': True,
# 			'add_ranksum_table': True
# 		}
# 	]
# }





# Retina models

CAN_LOAD_MODEL = True

# Mouse-only model
# MODEL_PATH = '/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/final_models/mouse_model.h5'
# LAYER_NAME = 'dense'
# FIT_DATA_KEYS = ['mouse_pos_val', 'mouse_neg_val']
# ACTIVATIONS_DIR = '/home/csestili/data/tacit_viz/retina_mouseonly/'

# # Multispecies model
MODEL_PATH = '/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/final_models/multispecies_model.h5'
LAYER_NAME = 'dense'
FIT_DATA_KEYS = ['mouse_pos_val', 'mouse_neg_val', 'human_pos_val', 'human_neg_val']
ACTIVATIONS_DIR = '/home/csestili/data/tacit_viz/retina_multispecies/'


DATA_MAPPING = {
	'mouse_pos_train': '/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/mouse4_v3/pos_TRAINING.fa',
	'mouse_pos_val': '/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/mouse4_v3/pos_VALIDATION.fa',
	'mouse_pos_test': '/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/mouse4_v3/pos_TEST.fa',
	'mouse_neg_train': '/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/mouse4_v3/neg_TRAINING.fa',
	'mouse_neg_val': '/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/mouse4_v3/neg_VALIDATION.fa',
	'mouse_neg_test': '/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/mouse4_v3/neg_TEST.fa',
	'human_pos_train': '/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/human4/pos_TRAINING.fa',
	'human_pos_val': '/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/human4/pos_VALIDATION.fa',
	'human_pos_test': '/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/human4/pos_TEST.fa',
	'human_neg_train': '/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/human4/neg_TRAINING.fa',
	'human_neg_val': '/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/human4/neg_VALIDATION.fa',
	'human_neg_test': '/projects/pfenninggroup/machineLearningForComputationalBiology/retina/models/human4/neg_TEST.fa',
	'dolphin': '/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/200_mammals/seqs/mm10orthologous_Tursiops_truncatus_500bp.fa',
	'rabbit': '/projects/pfenninggroup/machineLearningForComputationalBiology/retina/data/200_mammals/seqs/mm10orthologous_Oryctolagus_cuniculus_500bp.fa'
}

VISUALIZATION_MAPPING = {
	'fit_data': FIT_DATA_KEYS,
	'plot_data': [
		{
			'title': 'Mammals',
			'groups': [
				{'id': 5, 'name': 'Mouse +', 'sets': ['mouse_pos_train', 'mouse_pos_val', 'mouse_pos_test']},
				{'id': 4, 'name': 'Mouse -', 'sets': ['mouse_neg_train', 'mouse_neg_val', 'mouse_neg_test']},
				{'id': 3, 'name': 'Human +', 'sets': ['human_pos_train', 'human_pos_val', 'human_pos_test']},
				{'id': 2, 'name': 'Human -', 'sets': ['human_neg_train', 'human_neg_val', 'human_neg_test']},
				{'id': 1, 'name': 'Rabbit', 'sets': ['rabbit']},
				{'id': 0, 'name': 'Dolphin', 'sets': ['dolphin']}
			],
			'add_violinplot': True,
			'add_ranksum_table': True
		}
	]
}






def main(visualization_mapping=VISUALIZATION_MAPPING):

	# Compile list of datasets that need activations
	datasets = visualization_mapping['fit_data'] + [name for plot_spec in visualization_mapping['plot_data'] for group in plot_spec['groups'] for name in group['sets']]
	datasets = list(set(datasets))

	# Get activations for each dataset
	print("Getting activations...")
	activations = {}
	if not os.path.exists(ACTIVATIONS_DIR):
		os.makedirs(ACTIVATIONS_DIR)
	model = None
	for name in datasets:
		activations_path = os.path.join(ACTIVATIONS_DIR, f"{name}_activations.npy")
		if not os.path.exists(activations_path):
			if model is None:
				model = models.load_model(MODEL_PATH)
			path = DATA_MAPPING[name]

			tmp_data_save_path = None
			if not CAN_LOAD_MODEL:
				print(name)
				tmp_data_save_path = os.path.join(ACTIVATIONS_DIR, f"{name}_data.npy")
				if os.path.exists(tmp_data_save_path):
					continue
			activations[name] = models.get_activations(
				model, path, out_file=activations_path, layer_name=LAYER_NAME, tmp_data_save_path=tmp_data_save_path)
		else:
			activations[name] = np.load(activations_path)

	# Fit reducer
	reducer_path = os.path.join(ACTIVATIONS_DIR, f"{REDUCER_TYPE}_reducer.pkl")
	if not os.path.exists(reducer_path):
		fit_data = np.concatenate(
			[activations[name] for name in visualization_mapping['fit_data']],
			axis=0)
		fit_fn = visualization.umap_fit if REDUCER_TYPE == 'umap' else visualization.pca_fit
		reducer = fit_fn(fit_data, reducer_outfile=reducer_path)
	else:
		with open(reducer_path, 'rb') as f:
			reducer = pickle.load(f)
	# Variance explained: first 2 components as a fraction of all components
	variance_pc1 = reducer.explained_variance_[0] / np.sum(reducer.explained_variance_)
	variance_pc2 = reducer.explained_variance_[1] / np.sum(reducer.explained_variance_)

	# Transform and plot
	for plot_spec in visualization_mapping['plot_data']:
		# Get the relevant activations
		transform_data = [activations[name] for group in plot_spec['groups'] for name in group['sets']]
		transform_data = np.concatenate(transform_data, axis=0)
		transform_labels = [
			group['id'] for group in plot_spec['groups']
			for name in group['sets']
			for _ in range(len(activations[name]))]
		transform_labels = np.array(transform_labels)

		# Shuffle for visibility of all classes
		rng = np.random.default_rng()
		combined = np.concatenate((transform_data, np.expand_dims(transform_labels, axis=1)), axis=1)
		rng.shuffle(combined)
		transform_data, transform_labels = combined[:, :-1], combined[:, -1]

		# Transform
		print("Getting transformed samples...")
		title = plot_spec['title']
		transform_outfile = os.path.join(ACTIVATIONS_DIR, f"{title}_{REDUCER_TYPE}.npy")
		transform_label_outfile = os.path.join(ACTIVATIONS_DIR, f"{title}_{REDUCER_TYPE}_labels.npy")
		if not os.path.exists(transform_outfile):
			transformed = visualization.transform(reducer, transform_data, transform_outfile=transform_outfile)
			np.save(transform_label_outfile, transform_labels)
		else:
			transformed = np.load(transform_outfile)
			transform_labels = np.load(transform_label_outfile)

		# Visualize

		# Font settings
		plt.rc('font', size=3)



		plot_outfile = os.path.join(ACTIVATIONS_DIR, f"{title}_{REDUCER_TYPE}.{PLOT_EXTENSION}")
		label_mapping = {group['id']: group['name'] for group in plot_spec['groups']}
		add_histogram = plot_spec.get('add_histogram', False)
		add_violinplot = plot_spec.get('add_violinplot', False)
		add_ranksum_table = plot_spec.get('add_ranksum_table', True)
		fig, axs = visualization.scatter(transformed, plot_outfile=None, transform_labels=transform_labels,
			label_mapping=label_mapping, scatter_kwargs={"s": 0.7, "alpha": 0.7},
			add_histogram=add_histogram,
			add_violinplot=add_violinplot, add_ranksum_table=add_ranksum_table)
		fig.suptitle(f"{title}")
		ax_num = 0
		axs[ax_num].set_title("First 2 PCs")
		axs[ax_num].set_xlabel(f"PC 1 Variance explained: {variance_pc1:0.4}")
		axs[ax_num].set_ylabel(f"PC 2 Variance explained: {variance_pc2:0.4}")
		if add_histogram:
			ax_num += 1
			axs[ax_num].set_title("First PC")
			axs[ax_num].set_xlabel(f"PC 1")
			axs[ax_num].set_ylabel("Density")
		if add_violinplot:
			ax_num += 1
			axs[ax_num].set_title(f"First Principal Component. Variance explained: {variance_pc1:0.4}")
			axs[ax_num].set_ylabel("First Principal Component")

			# Put ticks inside plot
			axs[ax_num].tick_params(axis="y",direction="in", pad=-22, left="off", labelleft="on")
			axs[ax_num].tick_params(axis="x",direction="in", pad=-15, left="off", labelleft="on")

			# Set tick length
			axs[ax_num].tick_params(length=_cm_to_pt(0.058))


		if add_ranksum_table:
			ax_num += 1
			# The group with the highest ID is the one that gets compared to all the other groups
			group_names = sorted([(group['id'], group['name']) for group in plot_spec['groups']])
			ranksum_common_group_name = group_names[-1][1]
			axs[ax_num].set_title(f"Rank-sum statistics, {ranksum_common_group_name} vs. rest")



		# Re-scale figure
		# (width (in), height (in))
		size = fig.get_size_inches()
		sf = _cm_to_inch(FIG_WIDTH_CM) / size[0]
		fig.set_size_inches(sf * size[0], sf * size[1])


		# Save figure
		if PLOT_EXTENSION in ['png', 'both']:
			plot_outfile = os.path.join(ACTIVATIONS_DIR, f"{title}_{REDUCER_TYPE}.png")
			_savefig(plot_outfile, axs)
		if PLOT_EXTENSION in ['svg', 'both']:
			plot_outfile = os.path.join(ACTIVATIONS_DIR, f"{title}_{REDUCER_TYPE}.svg")
			_savefig(plot_outfile, axs)


def _savefig(plot_outfile, axs):
	if plot_outfile.endswith('svg'):
		# Clear points from scatterplot
		axs[0].clear()
	plt.savefig(plot_outfile, dpi=250)

def _cm_to_inch(len_cm):
	return len_cm * 0.393701

def _cm_to_pt(len_cm):
	return len_cm * 28.3465

if __name__ == '__main__':
	main()