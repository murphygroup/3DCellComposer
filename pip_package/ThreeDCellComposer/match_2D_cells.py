import numpy as np
from scipy.sparse import csr_matrix
from skimage.segmentation import find_boundaries
from datetime import datetime

"""
FUNCTIONS FOR MATCHING CELLS IN 2D SLICES
Author: Haoran Chen and Robert F. Murphy
Version: 1.3 February 14, 2025 R.F.Murphy, Ted Zhang
        Add option for setting JI_max
        Add logging
"""

def get_compartments_diff(arr1, arr2):
	a = set((tuple(i) for i in arr1))
	b = set((tuple(i) for i in arr2))
	diff = np.array(list(a - b))
	return diff


def get_matched_cells(current_cell_arr, new_cell_arr):
	a = set((tuple(i) for i in current_cell_arr))
	b = set((tuple(i) for i in new_cell_arr))
	JI = len(list(a & b)) / len(list(a | b))
	if JI != 0:
		return np.array(list(a)), np.array(list(b)), JI
	else:
		return False, False, False

def append_coord(rlabel_mask, indices, maxvalue):
	masked_imgs_coord = [[[], []] for i in range(maxvalue)]
	for i in range(0, len(rlabel_mask)):
		masked_imgs_coord[rlabel_mask[i]][0].append(indices[0][i])
		masked_imgs_coord[rlabel_mask[i]][1].append(indices[1][i])
	return masked_imgs_coord

def unravel_indices(labeled_mask, maxvalue):
	rlabel_mask = labeled_mask.reshape(-1)
	indices = np.arange(len(rlabel_mask))
	indices = np.unravel_index(indices, (labeled_mask.shape[0], labeled_mask.shape[1]))
	masked_imgs_coord = append_coord(rlabel_mask, indices, maxvalue)
	masked_imgs_coord = list(map(np.asarray, masked_imgs_coord))
	return masked_imgs_coord

def get_coordinates(mask):
	print("Getting cell coordinates...")
	cell_num = np.unique(mask)
	maxvalue = len(cell_num)
	channel_coords = unravel_indices(mask, maxvalue)
	return channel_coords

def compute_M(data):
	cols = np.arange(data.size)
	return csr_matrix((cols, (data.ravel(), cols)),
	                  shape=(data.max() + 1, data.size))

def get_indices_sparse(data):
	M = compute_M(data)
	return [np.unravel_index(row.data, data.shape) for row in M]

def list_remove(c_list, indexes):
	for index in sorted(indexes, reverse=True):
		del c_list[index]
	return c_list

def filter_cells(coords, mask):
	no_cells = []
	for i in range(len(coords)):
		if np.sum(mask[coords[i]]) == 0:
			no_cells.append(i)
	new_coords = list_remove(coords.copy(), no_cells)
	return new_coords

def get_indexed_mask(mask, boundary):
	boundary = boundary * 1
	boundary_loc = np.where(boundary == 1)
	boundary[boundary_loc] = mask[boundary_loc]
	return boundary

def get_boundary(mask):
	mask_boundary = find_boundaries(mask, mode='inner')
	mask_boundary_indexed = get_indexed_mask(mask, mask_boundary)
	return mask_boundary_indexed

def get_mask(cell_list):
	mask = np.zeros((img.shape))
	for cell_num in range(len(cell_list)):
		mask[tuple(cell_list[cell_num].T)] = cell_num
	return mask

def get_new_slice_mask(current_matched_index, new_matched, new_unmatched, max_cell_num, img_current_slice_shape):
	mask = np.zeros((img_current_slice_shape))
	for cell_num in range(len(new_matched)):
		mask[tuple(new_matched[cell_num].T)] = current_matched_index[cell_num]+1
	for cell_num in range(len(new_unmatched)):
		mask[tuple(new_unmatched[cell_num].T)] = max_cell_num + cell_num + 1
	return mask.astype(int)


def get_unmatched_list(matched_list, new_slice_cell_coords):
	total_num = len(new_slice_cell_coords)
	unmatched_list_index = []
	for i in range(total_num):
		if i not in matched_list:
			unmatched_list_index.append(i)
	unmatched_list = []
	for i in range(len(unmatched_list_index)):
		unmatched_list.append(new_slice_cell_coords[unmatched_list_index[i]])
	return unmatched_list, unmatched_list_index

def matching_cells_2D(img, JI_thre):
	#print(f"{datetime.now()} Starting matching_cells_2D with JI {JI_thre}...")
	new_img = [img[0]]
	#print(img.shape)
	for slice_num in range(1, img.shape[0]):
		#print("Matching", slice_num, 'th slice')
		img_current_slice = new_img[slice_num-1].astype(int)
		img_new_slice = img[slice_num].astype(int)
		current_slice_cell_coords = get_indices_sparse(img_current_slice)[1:]
		new_slice_cell_coords = get_indices_sparse(img_new_slice)[1:]
		
		current_slice_cell_coords = list(map(lambda x: np.array(x).T, current_slice_cell_coords))
		new_slice_cell_coords = list(map(lambda x: np.array(x).T, new_slice_cell_coords))
		
		#print(f"{datetime.now()} After get_indices_sparse")
		current_slice_cell_matched_list = []
		new_slice_cell_matched_list = []
		current_slice_cell_matched_index_list = []
		new_slice_cell_matched_index_list = []

		for i in range(len(current_slice_cell_coords)):
			if len(current_slice_cell_coords[i]) != 0:
				new_slice_search_num = np.unique(list(map(lambda x: img_new_slice[tuple(x)], current_slice_cell_coords[i])))
				best_JI = 0
				current_slice_cell_best = []
				for j in new_slice_search_num:
					if j != 0:
						if (j-1 not in new_slice_cell_matched_index_list) and (i not in current_slice_cell_matched_index_list):
							current_slice_cell, new_slice_cell, JI = get_matched_cells(current_slice_cell_coords[i], new_slice_cell_coords[j-1])
							if type(current_slice_cell) != bool:
								if JI > best_JI and JI > JI_thre:
									best_JI = JI
									current_slice_cell_best = current_slice_cell
									new_slice_cell_best = new_slice_cell
									i_ind = i
									j_ind = j-1
				
				if len(current_slice_cell_best) > 0:
					current_slice_cell_matched_list.append(current_slice_cell_best)
					new_slice_cell_matched_list.append(new_slice_cell_best)
					current_slice_cell_matched_index_list.append(i_ind)
					new_slice_cell_matched_index_list.append(j_ind)

		#print(f"{datetime.now()} After i,j loop")
		
		new_slice_cell_unmatched_list, new_slice_cell_unmatched_index_list = get_unmatched_list(new_slice_cell_matched_index_list, new_slice_cell_coords)
		new_slice_updated_mask = get_new_slice_mask(current_slice_cell_matched_index_list, new_slice_cell_matched_list, new_slice_cell_unmatched_list, len(np.unique(new_img[:slice_num])), img_current_slice.shape)
		new_img.append(new_slice_updated_mask)
		#print(f"{datetime.now()} After slice {slice_num}")
	new_img = np.stack(new_img, axis=0)
	#print(f"{datetime.now()} After matching")
	return new_img
