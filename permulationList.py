import pandas as pd
import argparse

'''
Used to generate peak list and filter prediction matrix for permulation 
Usage:

python permulationList.py -i longevity_perm_10k_computedP.csv -t 0.0005 -p 900000 \
-m $preds -o perm_1m_tmp

'''

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input', required=True, help="csv file with Exp_Pvalue column and a column with peak name")
	parser.add_argument('-t', '--threshold', required=True, help="threshold used to filter peaks by p-values")
	parser.add_argument('-p', '--permulation', required=False, help="value append to Missing_Trials column")
	parser.add_argument('-m', '--matrix', required=False, help="prediction matrix if wants to get filered, expect tsv, expect NO header")
	parser.add_argument('-n', '--colName', required=False, help="column name for p value")
	parser.add_argument('-o', '--output', required=True, help="output file name, e.g. perm10k")
	
	args = parser.parse_args()
	df = pd.read_csv(args.input, header='infer')
	if args.colName is not None:
		df = df.loc[df[args.colName] <= float(args.threshold)]
	else:
		df = df.loc[df['Exp_Pvalue'] <= float(args.threshold)] 

	if args.permulation is not None:
		df['Missing_Trials'] = args.permulation

	df.to_csv(args.output+"_peakList.csv", index=False)
	print("peak list output done")

	if args.matrix is not None:
		mat = pd.read_table(args.matrix, header=None)
		# Important: assumes peak name is in the first column
		filterMat = mat[mat.iloc[:, 0].isin(df.iloc[:, 0])]
		filterMat.to_csv(args.output+"_predictionMatrix.tsv", sep='\t', index=False, header=False)
		print("matrix filtering done")
