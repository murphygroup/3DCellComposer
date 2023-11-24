import os.path

import sys
from skimage.io import imread, imsave
import argparse
import tifffile
import numpy as np
from segmentation.segmentation_2D import deepcell_segmentation_2D
from segmentation.match_2D_cells import matching_cells_2D
from segmentation.match_3D_cells import matching_cells_3D
from segmentation.match_3D_nuclei import matching_nuclei_3D
from evaluation.single_method_eval_3D import seg_evaluation_3D
from visualization.coloring_3D import coloring_3D
from visualization.meshing_3D import meshing_3D
import pickle
import bz2

def get_channel_names(img_dir):
	with tifffile.TiffFile(img_dir) as tif:
		tif_tags = {}
		for tag in tif.pages[0].tags.values():
			name, value = tag.name, tag.value
			tif_tags[name] = value
	description = tif_tags['ImageDescription']
	name_list = list()
	for i in range(50):
		channel_num = "Channel:0:" + str(i)
		channel_anchor = description.find(channel_num)
		channel_str = description[channel_anchor:channel_anchor + 80]
		name_anchor = channel_str.find("Name")
		name_str = channel_str[name_anchor + 6:name_anchor + 20]
		channel_name = name_str[:name_str.find('"')]
		if len(channel_name) > 0:
			name_list.append(channel_name)
	return name_list


def get_channel_intensity(marker_list, names, img):
	channel_intensity = np.zeros((img.shape[0], img.shape[2], img.shape[3]))
	for marker in marker_list:
		channel_idx = names.index(marker)
		channel_intensity = channel_intensity + img[:, channel_idx, :, :]
	return channel_intensity

def parse_marker_list(arg):
	return arg.split(',')

def write_IMC_input_channels(img_dir, nucleus_channel_marker_list, cytoplasm_channel_marker_list, membrane_channel_marker_list):
	image = imread(img_dir)
	channel_names = get_channel_names(img_dir)
	
	# nucleus_channel_marker_list = ['Ir191']
	nucleus_channel = get_channel_intensity(nucleus_channel_marker_list, channel_names, image)
	
	# cytoplasm_channel_marker_list = ['In115', 'Y89', 'Tb159']
	cytoplasm_channel = get_channel_intensity(cytoplasm_channel_marker_list, channel_names, image)
	
	# membrane_channel_marker_list = ['La139', 'Pr141', 'Eu151', 'Gd160', 'Dy162']
	membrane_channel = get_channel_intensity(membrane_channel_marker_list, channel_names, image)
	
	imsave(f'{os.path.dirname(img_dir)}/nucleus.tif', nucleus_channel)
	imsave(f'{os.path.dirname(img_dir)}/cytoplasm.tif', nucleus_channel)
	imsave(f'{os.path.dirname(img_dir)}/membrane.tif', nucleus_channel)
	
	return nucleus_channel, cytoplasm_channel, membrane_channel, image

def main():
	parser = argparse.ArgumentParser(description="3DCellComposer Image Processor")
	parser.add_argument("image_path", help="Path to the multiplexed image")
	parser.add_argument("nucleus_channel_marker_list", type=parse_marker_list,
	                    help="A list of nuclear marker(s) in multiplexed image as input for segmentation")
	parser.add_argument("cytoplasm_channel_marker_list", type=parse_marker_list,
	                    help="A list of cytoplasmic marker(s) in multiplexed image as input for segmentation")
	parser.add_argument("membrane_channel_marker_list", type=parse_marker_list,
	                    help="A list of cell membrane marker(s) in multiplexed image as input for segmentation")
	
	args = parser.parse_args()
	# Process the image
	print("Generating input channels for segmentation...")
	nucleus_channel, cytoplasm_channel, membrane_channel, image = write_IMC_input_channels(args.image_path,
																						   args.nucleus_channel_marker_list,
																						   args.cytoplasm_channel_marker_list,
																						   args.membrane_channel_marker_list)
	
	print("Segmenting every 2D slices across three axes...")
	cell_mask_all_axes = {}
	nuclear_mask_all_axes = {}
	for axis in ['XY', 'XZ', 'YZ']:
		cell_mask_axis, nuclear_mask_axis = deepcell_segmentation_2D(nucleus_channel, membrane_channel, axis)
		cell_mask_all_axes[axis] = cell_mask_axis
		nuclear_mask_all_axes[axis] = nuclear_mask_axis
	
	print("Matching 2D cells in adjacent slices for each axis...")
	matched_2D_stack_all_JI = {}
	for JI in np.arange(0, 0.8, 0.1):
		matched_2D_stack_all_JI[JI] = {}
		for axis in ['XY', 'XZ', 'YZ']:
			matched_2D_stack_axis = matching_cells_2D(cell_mask_all_axes[axis], JI)
			matched_2D_stack_all_JI[JI][axis] = matched_2D_stack_axis
	
	print("Matching and repairing 3D cells...")
	matched_3D_all_JI = {}
	for JI in np.arange(0, 0.8, 0.1):
		matched_2D_stack_XY = matched_2D_stack_all_JI[JI]['XY']
		matched_2D_stack_XZ = matched_2D_stack_all_JI[JI]['XZ']
		matched_2D_stack_YZ = matched_2D_stack_all_JI[JI]['YZ']
		matched_3D_cell_mask = matching_cells_3D(matched_2D_stack_XY, matched_2D_stack_XZ, matched_2D_stack_YZ)
		matched_3D_all_JI[JI] = matched_3D_cell_mask
		
	print("Matching 3D nuclei...")
	final_matched_3D_cell_mask_JI = {}
	final_matched_3D_nuclear_mask_JI = {}
	for JI in np.arange(0, 0.8, 0.1):
		matched_3D_cell_mask = matched_3D_all_JI[JI]
		final_matched_3D_cell_mask, final_matched_3D_nuclear_mask = matching_nuclei_3D(matched_3D_cell_mask, nuclear_mask_all_axes['XY'])
		final_matched_3D_cell_mask_JI[JI] = final_matched_3D_cell_mask
		final_matched_3D_nuclear_mask_JI[JI] = final_matched_3D_nuclear_mask
	
	print("Evaluating 3D cell segmentation...")
	quality_score_JI = list()
	for JI in np.arange(0, 0.8, 0.1):
		final_matched_3D_cell_mask = final_matched_3D_cell_mask_JI[JI]
		final_matched_3D_nuclear_mask = final_matched_3D_nuclear_mask_JI[JI]
		quality_score = seg_evaluation_3D(final_matched_3D_cell_mask,
										  final_matched_3D_nuclear_mask,
										  image,
										  nucleus_channel,
										  cytoplasm_channel,
										  membrane_channel,
										  f'{os.path.dirname(os.path.dirname(args.image_path))}/evaluation/model')
		quality_score_JI.append(quality_score)
	best_quality_score = max(quality_score_JI)
	best_JI = quality_score_JI.index(best_quality_score)
	
	print(f"Quality Score of this 3D Cell Segmentation = {best_quality_score}")
	
	print("Generating surface meshes for visualization in Blender..")
	final_matched_3D_cell_mask = final_matched_3D_cell_mask_JI[best_JI]
	final_matched_3D_cell_mask_colored, number_of_colors = coloring_3D(final_matched_3D_cell_mask)
	meshing_3D(final_matched_3D_cell_mask, final_matched_3D_cell_mask_colored, number_of_colors)

	print("3D Segmentation and Evaluation Completed.")


if __name__ == "__main__":
	main()
