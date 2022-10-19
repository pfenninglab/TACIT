"""visualization.py: Tools for visualizing data."""
import pickle

import matplotlib.pyplot as plt
import numpy as np
import scipy
import sklearn
import umap

def get_umap(fit_data, fit_labels=None, transform_data=None, transform_labels=None, transform_outfile=None, plot_outfile=None, umap_kwargs=None, scatter_kwargs=None):
	"""Fit a UMAP reducer, and optionally transform and plot.

	Args:
		fit_data (np array), shape (num_samples, num_features): Data to fit the UMAP transform.
		fit_labels (np array), shape (num_samples,): If given, then use labels for
			supervised dimension reduction.
		transform_data (np array), shape (num_samples_tx, num_features): Data to transform
			after UMAP is fit. If not passed, then fit_data will be transformed.
		transform_labels (np array), shape (num_samples_tx,): If given, then use these labels
			to color the points in the visualization.
		transform_outfile (str): Path to save transformed data, .npy.
		plot_outfile (str): Path to save plot of transformed data, .jpg or .png.
		scatter_kwargs (dict): Keyword arguments to matplotlib scatter().
		umap_kwargs (dict): Keyword arguments to create the UMAP reducer.
	"""
	# TODO DRY

	# Fit UMAP transform
	umap_kwargs = umap_kwargs or {}
	reducer = umap.UMAP(**umap_kwargs)
	print("Fitting UMAP transform...")
	reducer.fit(fit_data, y=fit_labels)

	# Transform data
	if transform_data is None:
		transform_data = fit_data
	print("Transforming data...")
	transformed = reducer.transform(transform_data)
	if transform_outfile is not None:
		np.save(transform_outfile, transformed)

	# Plot transformed data
	if plot_outfile is not None:
		print("Plotting...")
		# Create default scatter_kwargs as an empty dict
		scatter_kwargs = scatter_kwargs or {}
		if transform_labels is not None:
			# Color the samples by their label
			scatter_kwargs['c'] = transform_labels
			scatter_kwargs['cmap'] = 'gist_rainbow'
		plt.scatter(transformed[:, 0], transformed[:, 1], **scatter_kwargs)
		if transform_labels is not None:
			# Add color legend
			num_labels = len(set(transform_labels))
			plt.colorbar(boundaries=np.arange(num_labels + 1) - 0.5).set_ticks(np.arange(num_labels))			
		plt.savefig(plot_outfile, dpi=300)

	return reducer, transformed

def umap_fit(fit_data, fit_labels=None, umap_kwargs=None, reducer_outfile=None):
	# Fit UMAP transform
	umap_kwargs = umap_kwargs or {}
	reducer = umap.UMAP(**umap_kwargs)
	print("Fitting UMAP transform...")
	reducer.fit(fit_data, y=fit_labels)

	# Save fit reducer object
	if reducer_outfile is not None:
		with open(reducer_outfile, 'wb') as f:
			pickle.dump(reducer, f)

	return reducer

def transform(reducer, transform_data, transform_outfile=None):
	print("Transforming data...")
	transformed = reducer.transform(transform_data)
	if transform_outfile is not None:
		np.save(transform_outfile, transformed)

	return transformed

def pca_fit(fit_data, pca_kwargs=None, reducer_outfile=None):
	# TODO DRY
	# Fit PCA transform
	pca_kwargs = pca_kwargs or {}
	reducer = sklearn.decomposition.PCA(**pca_kwargs)
	print("Fitting PCA transform...")
	reducer.fit(fit_data)
	print(f"PCA variance explained: {reducer.explained_variance_}")

	# Save fit reducer object
	if reducer_outfile is not None:
		with open(reducer_outfile, 'wb') as f:
			pickle.dump(reducer, f)

	return reducer

def scatter(points, plot_outfile=None, transform_labels=None, label_mapping=None, scatter_kwargs=None, add_histogram=False, add_violinplot=False, add_ranksum_table=False):
	print("Plotting...")
	# Create default scatter_kwargs as an empty dict
	scatter_kwargs = scatter_kwargs or {}
	if transform_labels is not None:
		unique_labels = np.unique(transform_labels)
		cmap = 'cool' if len(unique_labels) <= 2 else 'tab10'
		cmap = plt.cm.get_cmap(cmap)
		# Color the samples by their label
		scatter_kwargs['c'] = transform_labels
		scatter_kwargs['cmap'] = cmap

	plt.clf()
	nrows = sum([True, add_histogram, add_violinplot, add_ranksum_table])
	fig, axs = plt.subplots(nrows=nrows)
	ax_num = 0
	size = fig.get_size_inches()
	fig.set_size_inches(size[0], size[1] * nrows)

	# Scatter plot
	plot = axs[ax_num].scatter(points[:, 0], points[:, 1], **scatter_kwargs)
	if transform_labels is not None:
		# Convert numerical labels to string in the legend
		lines, labels = plot.legend_elements()
		if label_mapping is not None:
			# Assumes labels is sorted unique transform labels
			labels = [label_mapping[itm] for itm in sorted(set(transform_labels))]
		# Add color legend and format figure
		axs[ax_num].legend(lines, labels, loc='lower right', prop={'size': 5})
		axs[ax_num].tick_params(labelbottom=False, bottom=False, labelleft=False, left=False)

	# Histogram
	if add_histogram:
		ax_num += 1
		bins = np.linspace(np.min(points[:, 0]), np.max(points[:, 0]), num=32)
		for value in unique_labels:
			hist_points = points[np.where(transform_labels == value)]
			color = cmap((value - np.min(unique_labels))/(np.max(unique_labels) - np.min(unique_labels)))
			# Dataset is just PC 1 for each group
			axs[ax_num].hist(hist_points[:, 0], label=label_mapping[value], alpha=0.5, density=True, bins=bins, color=color)
			axs[ax_num].legend(loc='lower right', prop={'size': 5})
			axs[ax_num].tick_params(labelleft=False, left=False)
		# Share x-axis with scatterplot
		axs[ax_num].get_shared_x_axes().join(axs[0], axs[ax_num])

	# Violin plot
	if add_violinplot:
		ax_num += 1
		# Reverse sort so that we go from highest (positive) to lowest (negative)
		labels = sorted(unique_labels, reverse=True)
		# Dataset is just PC 1 for each group
		dataset = [points[np.where(transform_labels == value)] for value in labels]
		dataset = [group[:, 0] for group in dataset]
		violins = axs[ax_num].violinplot(dataset, showmeans=True)
		# We use the reversed label order for the xticks
		axs[ax_num].set_xticks(np.array(labels) + 1)
		# Don't need to reverse the labels, because they get mapped to the xticks in order
		axs[ax_num].set_xticklabels([label_mapping[itm] for itm in unique_labels], rotation=-70)
		# Color the violins
		for vp, value in zip(violins['bodies'], labels):
			color = cmap((value - np.min(unique_labels))/(np.max(unique_labels) - np.min(unique_labels)))
			vp.set_facecolor(color)
		# Set line thickness to 0.75pt
		for key in ['cmeans', 'cmaxes', 'cmins', 'cbars']:
			violins[key].set_linewidth(0.75)

	if add_ranksum_table:
		# Do rank-sum test between the group with the largest label, and every other group
		ax_num += 1
		# Reverse sort so that we go from highest (positive) to lowest (negative)
		labels = sorted(unique_labels, reverse=True)
		sample2 = points[np.where(transform_labels == labels[0])]
		table_data = [['Group', 'Statistic', 'P', 'N']]
		# Bonferroni correction, number of tests
		bonferroni_multiplier = len(labels) - 1
		for label in labels[1:]:
			sample1 = points[np.where(transform_labels == label)]
			ranksum_result = scipy.stats.ranksums(sample1[:, 0], sample2[:, 0])
			pvalue = ranksum_result.pvalue * bonferroni_multiplier
			table_data.append([label_mapping[label], ranksum_result.statistic, pvalue, len(sample1)])
		axs[ax_num].table(cellText=table_data, loc='center')
		axs[ax_num].axis('off')
		axs[ax_num].axis('tight')

	# Arrange plots so that they don't squeeze into each other
	plt.tight_layout(pad=5)
	if plot_outfile:
		plt.savefig(plot_outfile, dpi=250)
	return fig, axs
