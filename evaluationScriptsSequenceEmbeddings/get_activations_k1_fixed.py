import os

os.environ['THEANO_FLAGS'] = "mode=FAST_RUN,device=cuda0,force_device=True,floatX=float32"
from keras.models import load_model
import keras
import pandas as pd
import numpy as np

DATA_DIR = '/home/csestili/data/tacit_viz/rev2/'

# Multi-species motor cortext model
ARGS_MC = {
	'part': 'mc',
	'orth_species': ['Mouse', 'Macaque', 'Rat', 'Bat'],
	'model_path': '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/cortexEnhancerVsOtherSpeciesEnhMouseMacaqueRatBatPeakLessThan1kb_data_train.KerasModels/cortexEnhancerShort_multiSpecies_enhLooseNeg_500bp_conv5TotalTotalFiltNoInterPoolVeryHighMomL2SmallBatchPretrainBal_model.hdf5',
	'layer_name': 'activation_6',
}

def main(args):
	ortholog_mapping = pd.read_csv(os.path.join(DATA_DIR, args['part'], 'ortholog_mapping.csv'))
	model = None

	for idx, row in ortholog_mapping.iterrows():
		if row['Clade'] != 'Glires':
			continue
		print(idx)
		for orth_species in args['orth_species']:
			activations_path = os.path.join(DATA_DIR, args['part'], row['Species Name'], 'embeddings', orth_species + '_orth_embeddings.npy')
			if not os.path.exists(activations_path):
				if model is None:
					model = load_model(args['model_path'])
				_make_parent_dir(activations_path)
				npy_path = os.path.join(DATA_DIR, args['part'], row['Species Name'], 'seq_npy', orth_species + '_orthologs.npy')
				if os.path.exists(npy_path):
					get_activations(model, args['layer_name'], npy_path, activations_path)
				else:
					print("WARNING: input data " + npy_path + " does not exist, skipping")

def _make_parent_dir(path):
	"""Make parent directory of path"""
	if not os.path.exists(os.path.abspath(os.path.dirname(path))):
		os.makedirs(os.path.abspath(os.path.dirname(path)))

def get_activations(model, layer_name, data_path, out_path):
	out_layer = model.get_layer(layer_name)
	extractor = keras.models.Model(input=model.inputs, output=out_layer.output)
	sequences = np.load(data_path)
	sequences = sequences[:, np.newaxis, :, :].swapaxes(2, 3)
	res = extractor.predict(sequences, batch_size=512, verbose=1)
	np.save(out_path, res)

if __name__ == '__main__':
	main(ARGS_MC)
