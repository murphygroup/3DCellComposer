import numpy as np
import pickle
import bz2
import pandas as pd
import random
random.seed(12345)

def get_indices_pandas(data):
	d = data.ravel()
	f = lambda x: np.unravel_index(x.index, data.shape)
	return pd.Series(d).groupby(d).apply(f)

def find_neighbors(segmentation, i, j, k):
	# Get the label of the cell at the current voxel
	current_label = segmentation[i, j, k]
	
	# List of possible neighbors' coordinates
	neighbors_coords = [(i - 1, j, k), (i + 1, j, k),
						(i, j - 1, k), (i, j + 1, k),
						(i, j, k - 1), (i, j, k + 1)]
	
	# Retrieve the labels of the neighbors
	neighbors = set()
	for x, y, z in neighbors_coords:
		if 0 <= x < segmentation.shape[0] and 0 <= y < segmentation.shape[1] and 0 <= z < segmentation.shape[2]:
			if segmentation[x, y, z] != current_label:
				neighbors.add(segmentation[x, y, z])
	
	return neighbors


def color_cells(segmentation, cell_coords, cmap):
	
	# Dictionary to store color assignments
	color_assignments = {}
	
	# Define the colors, you can use more if needed
	cell_num = 1
	for label in cell_coords.index:
		# print(cell_num)
		cell_num += 1
		available_colors = set(cmap)
		
		
		# For each voxel in the segmentation, find neighbors and reduce available colors
		voxel_indices = np.vstack(cell_coords[label]).T
		for i, j, k in voxel_indices:
			neighbors = find_neighbors(segmentation, i, j, k)
			for neighbor in neighbors:
				if neighbor in color_assignments:
					available_colors.discard(color_assignments[neighbor])
		

		chosen_color = random.choice(list(available_colors))
		color_assignments[label] = chosen_color
		
	return color_assignments

def coloring(seg, coords, cols, cmap):
	seg_colored = np.zeros(seg.shape)
	for cell_idx in coords.index:
		seg_colored[coords[cell_idx]] = cmap.index(cols[cell_idx]) + 1
	return seg_colored

def get_large_color_list(number_of_colors):
	color_list = []
	max_value = 256
	for r in range(max_value):
		for g in range(max_value):
			for b in range(max_value):
				color_list.append((r, g, b))
				if len(color_list) == number_of_colors:
					return color_list
	return color_list



def coloring_3D(segmentation):
	cell_coords = get_indices_pandas(segmentation)[1:]
	color_codes = [1, 2, 3, 4]
	while True:
		try:
			segmentation_colored = color_cells(segmentation, cell_coords, color_codes)
			break  # If no error, exit the loop
		except Exception:
			print(f"Colors not enough. Adding a new color code and retrying.")
			color_codes.append(max(color_codes) + 1)
			

	return segmentation_colored, len(color_codes)
