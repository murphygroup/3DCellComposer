import numpy as np
from scipy.sparse import csr_matrix
import pandas as pd

"""
FUNCTIONS FOR MATCHING 3D CELLS
Author: Haoran Chen and Robert F. Murphy
Version: 1.3 February 14, 2025 R.F.Murphy, Ted Zhang
        remove rotations because already done in deepcell_segmentation_2D
"""

def compute_M(data):
	cols = np.arange(data.size)
	return csr_matrix((cols, (data.ravel(), cols)),
	                  shape=(data.max() + 1, data.size))

def get_indices_pandas(data):
	d = data.ravel()
	f = lambda x: np.unravel_index(x.index, data.shape)
	return pd.Series(d).groupby(d).apply(f)


def get_volumes_from_indices(series):
	return series.apply(lambda x: len(x[0]))
	
def matching_cells_3D(mask_XY, mask_XZ, mask_YZ):

	#print(mask_XY.shape,mask_XZ.shape,mask_YZ.shape)
	#this rotation is already done in "deepcell_segmentation_2D"
	#mask_XZ = np.rot90(mask_XZ, k=1, axes=(0, 2))
	#mask_YZ = np.rot90(mask_YZ, k=1, axes=(0, 1))

	#print(mask_XY.shape,mask_XZ.shape,mask_YZ.shape)

	X_max = np.max(mask_YZ) + 1
	Y_max = np.max(mask_XZ) + 1
	segmentation = np.zeros(mask_XY.shape, dtype=np.int64)

	
	for z in range(0, mask_XY.shape[0]):
		for x in range(0, mask_XY.shape[1]):
			for y in range(0, mask_XY.shape[2]):
				Z = mask_XY[z, x, y]
				Y = mask_XZ[z, x, y]
				X = mask_YZ[z, x, y]
				if Z == 0:
					index_1D = 0
				else:
					index_1D = Y + X * Y_max + Z * X_max * Y_max
				segmentation[z, x, y] = index_1D
	
	
	cell_coords = get_indices_pandas(segmentation)[1:]
	
	
	cell_volumes = get_volumes_from_indices(cell_coords)
	sorted_cells = sorted(cell_volumes.keys(), key=lambda x: cell_volumes[x], reverse=True)

	XY_coords = get_indices_pandas(mask_XY.astype(int))

	segmentation_XY_repaired = np.zeros(segmentation.shape, dtype=np.int64)
	segmentation_XY_repaired_binary = np.zeros(segmentation.shape)

	for cell_index in sorted_cells:
		current_coords = cell_coords[cell_index]
		if not np.any(segmentation_XY_repaired_binary[current_coords]) and len(np.unique(current_coords[0])) > 3:
			Z = int(cell_index // (X_max * Y_max))
			XY_coords_Z = XY_coords[Z]
			z_min = min(current_coords[0])
			z_max = max(current_coords[0])
			x_min = min(XY_coords_Z[1])
			x_max = max(XY_coords_Z[1])
			y_min = min(XY_coords_Z[2])
			y_max = max(XY_coords_Z[2])
			impute_range = np.where((XY_coords_Z[0] >= z_min) & (XY_coords_Z[0] <= z_max) & (XY_coords_Z[1] >= x_min) & (XY_coords_Z[1] <= x_max) & (XY_coords_Z[2] >= y_min) & (XY_coords_Z[2] <= y_max))
			XY_coords_Z = list(XY_coords_Z)
			XY_coords_Z[0] = XY_coords_Z[0][impute_range]
			XY_coords_Z[1] = XY_coords_Z[1][impute_range]
			XY_coords_Z[2] = XY_coords_Z[2][impute_range]
			XY_coords_Z = tuple(XY_coords_Z)
			segmentation_XY_repaired[XY_coords_Z] = cell_index
			segmentation_XY_repaired_binary[XY_coords_Z] = 1
		else:
			segmentation_XY_repaired_binary[current_coords] = 1
	segmentation_XY_repaired_output = np.zeros(segmentation_XY_repaired.shape, dtype=int)
	segmentation_XY_repaired_coords = get_indices_pandas(segmentation_XY_repaired)[1:]
	cell_idx = 1
	for coord_idx in segmentation_XY_repaired_coords.index:
		segmentation_XY_repaired_output[segmentation_XY_repaired_coords[coord_idx]] = cell_idx
		cell_idx += 1
	return segmentation_XY_repaired

