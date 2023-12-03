import os
import re
from typing import Any, Dict, List, Tuple
import numpy as np
import sys
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.stats import variation
from skimage.filters import threshold_mean, threshold_otsu
from skimage.morphology import area_closing, closing, disk
from skimage.segmentation import morphological_geodesic_active_contour as MorphGAC
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
from skimage.io import imread, imsave
import glob
import pickle
import bz2



def thresholding(img):
	threshold = threshold_mean(img.astype(np.int64))
	img_thre = img > threshold
	img_thre = img_thre * 1
	return img_thre


def fraction(img_bi, mask_bi):
	foreground_all = np.sum(img_bi)
	background_all = img_bi.shape[0] * img_bi.shape[1] * img_bi.shape[2] - foreground_all
	mask_all = np.sum(mask_bi)
	background = len(np.where(mask_bi - img_bi == 1)[0])
	foreground = np.sum(mask_bi * img_bi)
	if background_all == 0:
		background_fraction = 0
	else:
		background_fraction = background / background_all
	foreground_fraction = foreground / foreground_all
	mask_fraction = foreground / mask_all
	return foreground_fraction, background_fraction, mask_fraction


def foreground_separation(img_thre):

	contour_ref = img_thre.copy()

	img_thre = closing(img_thre, disk(1))

	img_thre = -img_thre + 1
	img_thre = closing(img_thre, disk(2))
	img_thre = -img_thre + 1


	img_thre = closing(img_thre, disk(10))

	img_thre = -img_thre + 1
	img_thre = closing(img_thre, disk(3))
	img_thre = -img_thre + 1


	img_thre = area_closing(img_thre, 5000, connectivity=2)


	contour_ref = contour_ref.astype(float)
	img_thre = img_thre.astype(float)
	img_binary = MorphGAC(
		-contour_ref + 1, 5, -img_thre + 1, smoothing=1, balloon=0.8, threshold=0.5
	)
	img_binary = area_closing(img_binary, 1000, connectivity=2)
	img_binary = -img_binary + 1

	return img_binary


def uniformity_CV(loc, channels):
	CV = []
	n = len(channels)
	for i in range(n):
		channel = channels[i]
		channel = channel / np.mean(channel)
		intensity = channel[tuple(loc.T)]
		CV.append(np.std(intensity))
	return np.average(CV)


def uniformity_fraction(loc, channels) -> float:
	n = len(channels)
	feature_matrix_pieces = []
	for i in range(n):
		channel = channels[i]
		ss = StandardScaler()
		z, x, y = channel.shape
		channel_z = ss.fit_transform(channel.reshape(z, x * y)).reshape(z, x, y)
		intensity = channel_z[tuple(loc.T)]
		feature_matrix_pieces.append(intensity)
	feature_matrix = np.vstack(feature_matrix_pieces)
	pca = PCA(n_components=1)
	model = pca.fit(feature_matrix.T)
	fraction = model.explained_variance_ratio_[0]
	return fraction


def foreground_uniformity(img_bi, mask, channels):
	foreground_loc = np.argwhere((img_bi - mask) == 1)
	CV = uniformity_CV(foreground_loc, channels)
	fraction = uniformity_fraction(foreground_loc, channels)
	return CV, fraction


def background_uniformity(img_bi, channels):
	background_loc = np.argwhere(img_bi == 0)
	CV = uniformity_CV(background_loc, channels)
	fraction = uniformity_fraction(background_loc, channels)
	return CV, fraction


def cell_uniformity_CV(feature_matrix):
	CV = []
	for i in range(feature_matrix.shape[1]):
		if np.sum(feature_matrix[:, i]) == 0:
			CV.append(np.nan)
		else:
			CV.append(variation(feature_matrix[:, i]))

	if np.sum(np.nan_to_num(CV)) == 0:
		return 0
	else:
		return np.nanmean(CV)


def cell_uniformity_fraction(feature_matrix):
	if np.sum(feature_matrix) == 0 or feature_matrix.shape[0] == 1:
		return 1
	else:
		pca = PCA(n_components=1)
		model = pca.fit(feature_matrix)
		fraction = model.explained_variance_ratio_[0]
		return fraction


def weighted_by_cluster(vector, labels):
	for i in range(len(vector)):
		vector[i] = vector[i] * len(np.where(labels == i)[0])
	weighted_average = np.sum(vector) / len(labels)
	return weighted_average


def cell_size_uniformity(mask):
	cell_coord = get_indices_pandas(mask)[1:]
	# cell_coord_num = len(cell_coord)
	cell_sizes = []
	for i in cell_coord.index:
		cell_size_current = len(cell_coord[i][0])
		if cell_size_current != 0:
			cell_sizes.append(cell_size_current)
	cell_sizes = np.expand_dims(np.array(cell_sizes), 1)
	cell_size_std = np.std(cell_sizes)
	cell_size_mean = np.mean(cell_sizes)
	cell_size_CV = cell_size_std / cell_size_mean
	return cell_size_CV, cell_sizes.T[0].tolist()


def cell_type(mask, channels):
	label_list = []
	n = len(channels)
	cell_coord = get_indices_pandas(mask)[1:]
	cell_coord_num = len(cell_coord)
	ss = StandardScaler()
	feature_matrix_z_pieces = []
	for i in range(n):
		channel = channels[i]
		z, x, y = channel.shape
		channel_z = ss.fit_transform(channel.reshape(z, x*y)).reshape(z, x, y)
		cell_intensity_z = []
		for j in cell_coord.index:
			cell_size_current = len(cell_coord[j][0])
			if cell_size_current != 0:
				single_cell_intensity_z = (
					np.sum(channel_z[tuple(cell_coord[j])]) / cell_size_current
				)
				cell_intensity_z.append(single_cell_intensity_z)
		feature_matrix_z_pieces.append(cell_intensity_z)

	feature_matrix_z = np.vstack(feature_matrix_z_pieces).T
	for c in range(1, 11):
		model = KMeans(n_clusters=c).fit(feature_matrix_z)
		label_list.append(model.labels_.astype(int))
	return label_list


def cell_uniformity(mask, channels, label_list):
	n = len(channels)
	cell_coord = get_indices_pandas(mask)[1:]
	cell_coord_num = len(cell_coord)
	ss = StandardScaler()
	feature_matrix_pieces = []
	feature_matrix_z_pieces = []
	for i in range(n):
		channel = channels[i]
		z, x, y = channel.shape
		channel_z = ss.fit_transform(channel.reshape(z, x*y)).reshape(z, x, y)
		cell_intensity = []
		cell_intensity_z = []
		for j in cell_coord.index:
			cell_size_current = len(cell_coord[j][0])
			if cell_size_current != 0:
				single_cell_intensity = np.sum(channel[tuple(cell_coord[j])]) / cell_size_current
				single_cell_intensity_z = (
					np.sum(channel_z[tuple(cell_coord[j])]) / cell_size_current
				)
				cell_intensity.append(single_cell_intensity)
				cell_intensity_z.append(single_cell_intensity_z)
		feature_matrix_pieces.append(cell_intensity)
		feature_matrix_z_pieces.append(cell_intensity_z)

	feature_matrix = np.vstack(feature_matrix_pieces).T
	feature_matrix_z = np.vstack(feature_matrix_z_pieces).T
	CV = []
	fraction = []
	silhouette = []

	for c in range(1, 11):
		labels = label_list[c - 1]
		CV_current = []
		fraction_current = []
		if c == 1:
			silhouette.append(1)
		else:
			silhouette.append(silhouette_score(feature_matrix_z, labels))
		for i in range(c):
			cluster_feature_matrix = feature_matrix[np.where(labels == i)[0], :]
			cluster_feature_matrix_z = feature_matrix_z[np.where(labels == i)[0], :]
			CV_current.append(cell_uniformity_CV(cluster_feature_matrix))
			fraction_current.append(cell_uniformity_fraction(cluster_feature_matrix_z))
		CV.append(weighted_by_cluster(CV_current, labels))
		fraction.append(weighted_by_cluster(fraction_current, labels))
	return CV, fraction, silhouette


def compute_M(data):
	cols = np.arange(data.size)
	return csr_matrix((cols, (data.ravel(), cols)), shape=(np.int64(data.max() + 1), data.size))


def get_indices_sparse(data):
	data = data.astype(np.uint64)
	M = compute_M(data)
	return [np.unravel_index(row.data, data.shape) for row in M]

def get_indices_pandas(data):
	d = data.ravel()
	f = lambda x: np.unravel_index(x.index, data.shape)
	return pd.Series(d).groupby(d).apply(f)


def get_indexed_mask(mask, boundary):
	boundary = boundary * 1
	boundary_loc = np.where(boundary == 1)
	boundary[boundary_loc] = mask[boundary_loc]
	return boundary


def flatten_dict(input_dict):
	local_list = []
	for key, value in input_dict.items():
		if type(value) == dict:
			local_list.extend(flatten_dict(value))
		else:
			local_list.append(value)
	return local_list



def get_quality_score(features, model):
	ss = model[0]
	pca = model[1]
	features_scaled = ss.transform(features)
	score = (
		pca.transform(features_scaled)[0, 0] * pca.explained_variance_ratio_[0]
		+ pca.transform(features_scaled)[0, 1] * pca.explained_variance_ratio_[1]
	)
	return score





def seg_evaluation_3D(cell_matched_mask,
                      nuclear_matched_mask,
                      nucleus,
                      cytoplasm,
                      membrane,
                      img_channels,
                      voxel_size,
                      pca_dir):

	
	cell_outside_nucleus_mask = cell_matched_mask - nuclear_matched_mask
	

	metric_mask = np.expand_dims(cell_matched_mask, 0)
	metric_mask = np.vstack((metric_mask, np.expand_dims(nuclear_matched_mask, 0)))
	metric_mask = np.vstack((metric_mask, np.expand_dims(cell_outside_nucleus_mask, 0)))

	if not os.path.exists(f'{pca_dir}/img_binary.pkl'):
		img_input_channel = nucleus + cytoplasm + membrane
		img_thresholded = np.stack([thresholding(slice_2d) for slice_2d in img_input_channel], axis=0)
	
		img_thresholded = np.sign(img_thresholded)
	
		img_binary_pieces = []
		for z in range(img_thresholded.shape[0]):
			img_binary_pieces.append(foreground_separation(img_thresholded[z]))
		img_binary = np.stack(img_binary_pieces, axis=0)
		img_binary = np.sign(img_binary)
		with bz2.BZ2File(f'{pca_dir}/img_binary.pkl', 'wb') as f:
			pickle.dump(img_binary, f)
	else:
		with bz2.BZ2File(f'{pca_dir}/img_binary.pkl', 'rb') as f:
			img_binary = pickle.load(f)
		
	# set mask channel names
	channel_names = [
		"Matched Cell",
		"Nucleus (including nuclear membrane)",
		"Cell Not Including Nucleus (cell membrane plus cytoplasm)",
	]
	metrics = {}

	# 3D IMC has dim order ZCXY
	img_channels = np.transpose(img_channels, (1, 0, 2, 3))

	if len(np.unique(cell_matched_mask))-1 > 10:
		print(f'- {len(np.unique(cell_matched_mask))} 3D cells segmented')
		for channel in range(metric_mask.shape[0]):
			current_mask = metric_mask[channel]
			mask_binary = np.sign(current_mask)
			metrics[channel_names[channel]] = {}
			if channel_names[channel] == "Matched Cell":
				voxel_volume = float(voxel_size[0]) * float(voxel_size[1]) * float(voxel_size[2])
				pixel_num = mask_binary.shape[0] * mask_binary.shape[1] * mask_binary.shape[2]
				micron_num = voxel_volume * pixel_num
	
				# TODO: match 3D cell and nuclei and calculate the fraction of match, assume cell and nuclei are matched for now
	
				# calculate number of cell per 100 cubic micron
				cell_num = len(np.unique(current_mask)) - 1
	
				cell_num_normalized = cell_num / micron_num * 100
	
				# calculate the standard deviation of cell size
				cell_size_CV, cell_sizes_voxels = cell_size_uniformity(current_mask)
				
				cell_sizes_microns = [size * voxel_volume for size in cell_sizes_voxels]
				weighted_avg_microns = sum(size * size for size in cell_sizes_microns) / sum(cell_sizes_microns)
				
				# get coverage metrics
				foreground_fraction, background_fraction, mask_foreground_fraction = fraction(
					img_binary, mask_binary
				)
				
				foreground_CV, foreground_PCA = foreground_uniformity(
					img_binary, mask_binary, img_channels
				)
				metrics[channel_names[channel]][
					"NumberOfCellsPer100CubicMicrons"
				] = cell_num_normalized
				metrics[channel_names[channel]][
					"FractionOfForegroundOccupiedByCells"
				] = foreground_fraction
				metrics[channel_names[channel]]["1-FractionOfBackgroundOccupiedByCells"] = (
						1 - background_fraction
				)
				metrics[channel_names[channel]][
					"FractionOfCellMaskInForeground"
				] = mask_foreground_fraction
				metrics[channel_names[channel]]["1/(CVOfCellSize+1)"] = 1 / (
						cell_size_CV + 1
				)
				
				metrics[channel_names[channel]][
					"WeightedAvgCellSizeinCubicMicrons"
				] = weighted_avg_microns
				
				metrics[channel_names[channel]]["1/(AvgCVForegroundOutsideCells+1)"] = 1 / (
						foreground_CV + 1
				)
				metrics[channel_names[channel]][
					"FractionOfFirstPCForegroundOutsideCells"
				] = foreground_PCA

				cell_type_labels = cell_type(current_mask, img_channels)
				
			else:
				# get cell uniformity
				cell_CV, cell_fraction, cell_silhouette = cell_uniformity(
					current_mask, img_channels, cell_type_labels
				)
				avg_cell_CV = np.average(cell_CV[0])
				avg_cell_fraction = np.average(cell_fraction[0])
				avg_cell_silhouette = np.average(cell_silhouette)
				metrics[channel_names[channel]][
					"1/(AvgOfWeightedAvgCVMeanCellIntensitiesOver1~10NumberOfClusters+1)"
				] = 1 / (avg_cell_CV + 1)
				metrics[channel_names[channel]][
					"AvgOfWeightedAvgFractionOfFirstPCMeanCellIntensitiesOver1~10NumberOfClusters"
				] = avg_cell_fraction
				metrics[channel_names[channel]][
					"AvgSilhouetteOver1~10NumberOfClusters"
				] = avg_cell_silhouette
	

		
		metrics_flat = np.expand_dims(flatten_dict(metrics), 0)
		ss, pca_model = pickle.load(open(f'{pca_dir}/pca_tissue.pkl', 'rb'))
		metrics_flat_z = ss.transform(metrics_flat)
		metrics_pca = pca_model.transform(metrics_flat_z)
		weighted_score = np.exp(sum(metrics_pca[0, i] * pca_model.explained_variance_ratio_[i] for i in range(2)))
		
		return weighted_score
	

