import numpy as np
from skimage.segmentation import find_boundaries
from scipy.sparse import csr_matrix
import pandas as pd

def match_3D_slice(mask_3D, cell_mask_3D):
	for z in range(mask_3D.shape[0]):
		for y in range(mask_3D.shape[1]):
			for x in range(mask_3D.shape[2]):
				if mask_3D[z, y, x] != 0:
					mask_3D[z, y, x] = cell_mask_3D[z, y, x]
	return mask_3D


def get_matched_cells(cell_arr, cell_membrane_arr, nuclear_arr, mismatch_repair):
	a = set((tuple(i) for i in cell_arr))
	b = set((tuple(i) for i in cell_membrane_arr))
	c = set((tuple(i) for i in nuclear_arr))
	d = a - b
	mismatch_pixel_num = len(list(c - d))
	# print(mismatch_pixel_num)
	mismatch_fraction = len(list(c - d)) / len(list(c))
	if not mismatch_repair:
		if mismatch_pixel_num == 0:
			return np.array(list(a)), np.array(list(c)), 0
		else:
			return False, False, False
	else:
		if mismatch_pixel_num < len(c):
			return np.array(list(a)), np.array(list(d & c)), mismatch_fraction
		else:
			return False, False, False


def compute_M(data):
	cols = np.arange(data.size)
	return csr_matrix((cols, (data.ravel(), cols)),
	                  shape=(data.max() + 1, data.size))

def get_indices_pandas(data):
	M = compute_M(data)
	return [np.unravel_index(row.data, data.shape) for row in M]


def get_indexed_mask(mask, boundary):
	boundary = boundary * 1
	boundary_loc = np.where(boundary == 1)
	boundary[boundary_loc] = mask[boundary_loc]
	return boundary

def get_boundary(mask):
	mask_boundary = find_boundaries(mask, mode='inner')
	mask_boundary_indexed = get_indexed_mask(mask, mask_boundary)
	return mask_boundary_indexed

def get_mask(cell_list, mask_shape):
	mask = np.zeros((mask_shape))
	for cell_num in range(len(cell_list)):
		mask[tuple(cell_list[cell_num].T)] = cell_num
	return mask

def get_indices_pandas(data):
	d = data.ravel()
	f = lambda x: np.unravel_index(x.index, data.shape)
	return pd.Series(d).groupby(d).apply(f)

	
def match_repair_cell_nucleus(nuclear_slice, cell_slice):
	cell_membrane_mask = get_boundary(cell_slice)
	cell_coords = get_indices_pandas(cell_slice)[1:]
	nucleus_coords = get_indices_pandas(nuclear_slice)[1:]
	cell_membrane_coords = get_indices_pandas(cell_membrane_mask)[1:]

	cell_matched_index_list = []
	nucleus_matched_index_list = []
	cell_matched_list = []
	nucleus_matched_list = []
	repaired_num = 0
	mismatch_repair = True
	for i in cell_coords.index:
		if len(cell_coords[i]) != 0 and np.unique(nuclear_slice)[0] != 1:
			current_cell_coords = np.array(cell_coords[i]).T
			current_cell_membrane_coords = np.array(cell_membrane_coords[i]).T
			nuclear_search_num = np.unique(list(map(lambda x: nuclear_slice[tuple(x)], current_cell_coords))).astype(np.int64)
			best_mismatch_fraction = 1
			whole_cell_best = []
			for j in nuclear_search_num:
				if j != 0:
					if (j not in nucleus_matched_index_list) and (i not in cell_matched_index_list):
						current_nucleus_coords = np.array(nucleus_coords[j]).T
						whole_cell, nucleus, mismatch_fraction = get_matched_cells(current_cell_coords, current_cell_membrane_coords, current_nucleus_coords, mismatch_repair=mismatch_repair)
						if type(whole_cell) != bool:
							if mismatch_fraction < best_mismatch_fraction:
								best_mismatch_fraction = mismatch_fraction
								whole_cell_best = whole_cell
								nucleus_best = nucleus
								j_ind = j
			if best_mismatch_fraction < 1 and best_mismatch_fraction > 0:
				repaired_num += 1
			
			if len(whole_cell_best) > 0:
				nucleus_matched_list.append(nucleus_best)
				nucleus_matched_index_list.append(j_ind)
	
	nuclear_matched_mask = get_mask(nucleus_matched_list, nuclear_slice.shape)

	return nuclear_matched_mask

def trim_nuclei_z_slice(cell_mask, nuclear_mask):
	matched_cell_mask = np.zeros(cell_mask.shape)
	matched_nuclear_mask = np.zeros(cell_mask.shape)
	cell_coords = get_indices_pandas(cell_mask)[1:]
	nuclear_coords = get_indices_pandas(nuclear_mask)[1:]
	for cell_idx in cell_coords.index:
		current_cell_coords = cell_coords[cell_idx]
		nuclear_all_zeros = not np.any(nuclear_mask[current_cell_coords])
		if not nuclear_all_zeros:
			current_nucleus_coords = nuclear_coords[cell_idx]
			current_nucleus_coords = np.array(current_nucleus_coords).T
			cell_start_slice_idx = np.min(current_cell_coords[0])
			cell_end_slice_idx = np.max(current_cell_coords[0])
			nuclear_start_slice_array_idx = np.where(current_nucleus_coords[:, 0] == cell_start_slice_idx)
			current_nucleus_coords = np.delete(current_nucleus_coords, nuclear_start_slice_array_idx, axis=0)
			nuclear_end_slice_array_idx = np.where(current_nucleus_coords[:, 0] == cell_end_slice_idx)
			current_nucleus_coords = np.delete(current_nucleus_coords, nuclear_end_slice_array_idx, axis=0)
			if current_nucleus_coords.shape[0] > 0:
				current_nucleus_coords = tuple(current_nucleus_coords.T)
				matched_nuclear_mask[current_nucleus_coords] = cell_idx
				matched_cell_mask[current_cell_coords] = cell_idx

			
	return matched_cell_mask, matched_nuclear_mask
	
	

def matching_nuclei_3D(cell_mask_3D, nuclear_slice):
	

	repaired_nuclear_slices = []
	for slice_idx in range(nuclear_slice.shape[0]):
		#print(slice_idx)
		repaired_nuclear_slice = match_repair_cell_nucleus(nuclear_slice[slice_idx], cell_mask_3D[slice_idx])
		repaired_nuclear_slices.append(repaired_nuclear_slice)
	
	repaired_nuclear_slice = np.stack(repaired_nuclear_slices)
	matched_nuclear_3D = match_3D_slice(repaired_nuclear_slice, cell_mask_3D)
	
	matched_cell_3D, matched_nuclear_3D = trim_nuclei_z_slice(cell_mask_3D, matched_nuclear_3D)
	return matched_cell_3D, matched_nuclear_3D
	
