import numpy as np
import bz2
import pickle
from pathlib import Path
from skimage import measure
from skimage.io import imsave
import pandas as pd
import colorsys

"""
WRAPPER TO GENERATE BLENDER FILES FOR VISUALIZATION
Author: Haoran Chen
Version: 1.1 December 14, 2023 Haoran Chen
        Fix output dir and update color map to 0-1
"""


def write_to_mtl(materials, filename):
	with open(filename, 'w') as f:
		for material_id, color in materials.items():

			f.write(f'newmtl Color_{material_id}\n')
			f.write(f'Kd {color[0]} {color[1]} {color[2]}\n')  # Diffuse color


def write_to_obj(verts, faces, groups, colors, filename):
	with open(filename, 'w') as f:
		f.write(
			f'mtllib ./data/cell_mesh.mtl\n')
		
		# Write vertices
		for v in verts:
			f.write('v {0} {1} {2}\n'.format(v[0], v[1], v[2]))
		
		# Keep track of the last group number to identify when to change groups
		last_group = None
		
		# Write faces
		for face_idx in range(len(faces)):
			# Take the value from the first vertex in the face as the group number
			face = faces[face_idx]
			group = groups[face_idx][0]
			color = colors[face_idx][0]
			
			# If this is a new group, write a new group tag
			if group != last_group:
				f.write(f'g Cell_{group}\n')
				last_group = group
				# print(group)
				f.write(f'usemtl Color_{color}\n')
			
			# Write the face
			f.write('f {0} {1} {2}\n'.format(face[0] + 1, face[1] + 1, face[2] + 1))




def get_indices_pandas(data):
	d = data.ravel()
	f = lambda x: np.unravel_index(x.index, data.shape)
	return pd.Series(d).groupby(d).apply(f)

# Extract the 2D contours (meshes) from the slices
def get_2D_mesh(slice_data, level=1.9):
	return measure.find_contours(slice_data, level=level)

# Convert 2D contours to 3D
def convert_2D_contour_to_3D(contour, z_value):
	return np.hstack((np.full((contour.shape[0], 1), z_value), contour))

# Triangulate the 2D contours to create triangles from them
def triangulate_2D_contour(contour):
	triangles = []
	for i in range(1, len(contour) - 1):
		triangles.append([0, i, i + 1])
	return triangles

import colorsys

def generate_color_map(number_of_colors):
	color_map = {}
	step = 360 / number_of_colors  # Divide the color wheel into equal parts

	for i in range(1, number_of_colors + 1):
		hue = (i - 1) * step
		# Convert HSL to RGB. Saturation and Lightness are set to 0.5 (50%) for vivid colors
		rgb_color = colorsys.hls_to_rgb(hue / 360, 0.5, 0.5)
		color_map[i] = rgb_color
	
	return color_map



def meshing_3D(mask, mask_colored, num_of_col, output_path: Path):

	cell_coords = get_indices_pandas(mask)[1:]
	
	
	
	all_verts = []
	all_faces = []
	all_values = []
	all_groups = []
	all_colors = []
	offset = 0  # To keep track of the index offset for faces when combining multiple cells
	
	for cell_index in cell_coords.index:
		current_coords = cell_coords[cell_index]
		current_mask = np.zeros(mask.shape)
		current_mask[current_coords] = 2
		
		# 3D mesh
		current_color = mask_colored[cell_index]
		verts, faces, normals, values = measure.marching_cubes(current_mask, level=1.999)
		faces += offset
		offset += len(verts)

		# 2D mesh for start and end slice to close the hole
		if len(np.unique(current_coords[0])) >= 2:
			z_start = np.min(current_coords[0])
			if np.sum(current_mask[z_start] != 0) > np.sum(current_mask[z_start+1] != 0):
				start_slice_mesh = get_2D_mesh(current_mask[z_start])
				start_2D_contours = [convert_2D_contour_to_3D(contour, z_start) for contour in start_slice_mesh]
				start_triangles = [triangulate_2D_contour(contour) for contour in start_2D_contours]
				start_2D_contours = np.vstack(start_2D_contours)
				start_triangles = np.vstack(start_triangles)
				verts = np.vstack([verts, start_2D_contours])
				start_triangles += offset
				offset += len(start_2D_contours)
				faces = np.vstack([faces, start_triangles])
			
			z_end = np.max(current_coords[0])
			if np.sum(current_mask[z_end] != 0) > np.sum(current_mask[z_end-1] != 0):
				end_slice_mesh = get_2D_mesh(current_mask[z_end])
				end_2D_contours = [convert_2D_contour_to_3D(contour, z_end) for contour in end_slice_mesh]
				end_triangles = [triangulate_2D_contour(contour) for contour in end_2D_contours]
				end_2D_contours = np.vstack(end_2D_contours)
				end_triangles = np.vstack(end_triangles)
				verts = np.vstack([verts, end_2D_contours])
				end_triangles += offset
				offset += len(end_2D_contours)
				faces = np.vstack([faces, end_triangles])
		
		# Append the new vertices, faces, and values to the main list
		all_verts.extend(verts)
		all_faces.extend(faces)
		all_values.extend(values)
		all_groups.extend(np.repeat(cell_index, len(faces)))
		all_colors.extend(np.repeat(current_color, len(faces)))

	all_verts = np.vstack(all_verts)
	all_faces = np.vstack(all_faces)
	all_groups = np.vstack(all_groups)
	all_colors = np.vstack(all_colors)
	
	color_map = generate_color_map(num_of_col)
	write_to_mtl(color_map, f'{output_path}/cell_mesh.mtl')
	write_to_obj(all_verts, all_faces, all_groups, all_colors, f'{output_path}/cell_mesh.obj')
