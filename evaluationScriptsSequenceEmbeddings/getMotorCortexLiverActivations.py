import os

os.environ['THEANO_FLAGS'] = "mode=FAST_RUN,device=cuda0,force_device=True,floatX=float32"
from keras.models import load_model
import keras
import numpy as np

## Motor cortex multispecies
#MODEL_PATH = '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/cortexEnhancerVsOtherSpeciesEnhMouseMacaqueRatBatPeakLessThan1kb_data_train.KerasModels/cortexEnhancerShort_multiSpecies_enhLooseNeg_500bp_conv5TotalTotalFiltNoInterPoolVeryHighMomL2SmallBatchPretrainBal_model.hdf5'
#LAYER_NAME = 'activation_6'
#FIT_DATA_KEYS = ['mouse_pos_val', 'mouse_neg_val',
#	'mac_pos_val', 'rat_pos_val', 'bat_pos_val',
#	'mac_neg_val', 'rat_neg_val', 'bat_neg_val']
#ACTIVATIONS_DIR = '/home/csestili/data/tacit_viz/mc_multispecies/'

## Motor cortex mouse-only
#MODEL_PATH = '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/cortexEnhancerVsOtherSpeciesEnhLoosePeakLessThan1kb500bp_data_train.KerasModels/hybrid_cortexEnhancerShort_enhLooseNeg_500bp_conv5MuchMuchMoreFiltNoInterPoolTwoDenseLargeDenseVeryHighMomLowDropL2SmallBatchPretrainBal_model.hdf5'
## flatten -> dense -> activation_6
#LAYER_NAME = 'activation_6'
#FIT_DATA_KEYS = ['mouse_pos_val', 'mouse_neg_val']
#ACTIVATIONS_DIR = '/home/csestili/data/tacit_viz/mc_mouseonly/'

#VISUALIZATION_MAPPING = {
#	'fit_data': FIT_DATA_KEYS,
#	'plot_data': [
#		{
#			'title': 'Mammals (order 2)',
#			'groups': [
#				{'id': 9, 'name': 'Mouse +', 'sets': ['mouse_pos']},
#				{'id': 8, 'name': 'Mouse -', 'sets': ['mouse_neg']},
#				{'id': 7, 'name': 'Macaque +', 'sets': ['mac_pos_train', 'mac_pos_val', 'mac_pos_test']},
#				{'id': 6, 'name': 'Macaque -', 'sets': ['mac_neg_train', 'mac_neg_val', 'mac_neg_test']},
#				{'id': 5, 'name': 'Rat +', 'sets': ['rat_pos_train', 'rat_pos_val', 'rat_pos_test']},
#				{'id': 4, 'name': 'Rat -', 'sets': ['rat_neg_train', 'rat_neg_val', 'rat_neg_test']},
#				{'id': 3, 'name': 'Bat +', 'sets': ['bat_pos_train', 'bat_pos_val', 'bat_pos_test']},
#				{'id': 2, 'name': 'Bat -', 'sets': ['bat_neg_train', 'bat_neg_val', 'bat_neg_test']},
#				{'id': 1, 'name': 'Rabbit', 'sets': ['rabbit']},
#				{'id': 0, 'name': 'Dolphin', 'sets': ['dolphin']}
#			],
#			'add_violinplot': True,
#			'add_ranksum_table': True
#		}
#	]
#}



# liver model
MODEL_PATH = '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/liverEnhancerVsOtherSpeciesEnhMouseMacaqueRatCowPigPeakLessThan1kb_data_train.KerasModels/liver_mouseMacaqueRatCowPig_otherSpeciesEnhLooseNeg_500bp_conv5AllFiltNoInterPoolVeryHighMomL2SmallBatchPretrainBal_model.hdf5'
LAYER_NAME = 'activation_6'
FIT_DATA_KEYS = ['mouse_pos_val', 'mouse_neg_val',
	'mac_pos_val', 'rat_pos_val', 'cow_pos_val', 'pig_pos_val',
	'mac_neg_val', 'rat_neg_val', 'cow_neg_val', 'pig_neg_val']
ACTIVATIONS_DIR = '/home/csestili/data/tacit_viz/liver_multispecies/'
VISUALIZATION_MAPPING = {
	'fit_data': FIT_DATA_KEYS,
	'plot_data': [
		{
			'title': 'Mammals',
			'groups': [
				{'id': 11, 'name': 'Mouse +', 'sets': ['mouse_pos']},
				{'id': 10, 'name': 'Mouse -', 'sets': ['mouse_neg']},
				{'id': 9, 'name': 'Macaque +', 'sets': ['mac_pos_train', 'mac_pos_val', 'mac_pos_test']},
				{'id': 8, 'name': 'Macaque -', 'sets': ['mac_neg_train', 'mac_neg_val', 'mac_neg_test']},
				{'id': 7, 'name': 'Rat +', 'sets': ['rat_pos_train', 'rat_pos_val', 'rat_pos_test']},
				{'id': 6, 'name': 'Rat -', 'sets': ['rat_neg_train', 'rat_neg_val', 'rat_neg_test']},
				{'id': 5, 'name': 'Cow +', 'sets': ['cow_pos_train', 'cow_pos_val', 'cow_pos_test']},
				{'id': 4, 'name': 'Cow -', 'sets': ['cow_neg_train', 'cow_neg_val', 'cow_neg_test']},
				{'id': 3, 'name': 'Pig +', 'sets': ['pig_pos_train', 'pig_pos_val', 'pig_pos_test']},
				{'id': 2, 'name': 'Pig -', 'sets': ['pig_neg_train', 'pig_neg_val', 'pig_neg_test']},
				{'id': 1, 'name': 'Rabbit', 'sets': ['rabbit']},
				{'id': 0, 'name': 'Dolphin', 'sets': ['dolphin']}
			],
			'add_violinplot': True,
			'add_ranksum_table': True
		}
	]
}



def main(visualization_mapping=VISUALIZATION_MAPPING):
	# Just get activations

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
				model = load_model(MODEL_PATH)
			print(name)
			data_path = os.path.join(ACTIVATIONS_DIR, f"{name}_data.npy")
			get_activations(model, LAYER_NAME, data_path, activations_path)

def get_activations(model, layer_name, data_path, out_path):
	out_layer = model.get_layer(layer_name)
	extractor = keras.models.Model(input=model.inputs, output=out_layer.output)
	sequences = np.load(data_path)
	sequences = sequences[:, np.newaxis, :, :].swapaxes(2, 3)
	res = extractor.predict(sequences, batch_size=512, verbose=1)
	np.save(out_path, res)

if __name__ == '__main__':
	main()
