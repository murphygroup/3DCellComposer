import os.path
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
import json

from ThreeDCellComposer.preprocessing import *
from ThreeDCellComposer.deepcell_only import deepcell_segmentation_2D
#from ThreeDCellComposer.custom_segmentation_wrapper import custom_segmentation_wrapper
#from .installation.install_all_methods import install_segmentation_methods
#from segmentation_2D.all_methods_segmentation import *
from ThreeDCellComposer.match_2D_cells import matching_cells_2D
from ThreeDCellComposer.match_3D_cells import matching_cells_3D
from ThreeDCellComposer.match_3D_nuclei import matching_nuclei_3D
from ThreeDCellComposer.single_method_eval_3D import seg_evaluation_3D
from ThreeDCellComposer.coloring_3D import coloring_3D
from ThreeDCellComposer.meshing_3D import meshing_3D

"""
ThreeDCellComposer MAIN PROGRAM
Author: Haoran Chen and Robert F. Murphy
Version: 1.1 December 14, 2023 R.F.Murphy, Haoran Chen
        Fixes in
        deepcell_only.py
        single_method_eval_3D.py
        meshing_3D.py
         1.2 February 8, 2024 R.F.Murphy
        created PyPI package
         1.5.1 May 27, 2025 R.F.Murphy
        added downsample support
"""

def parse_marker_list(arg):
	return arg.split(',')


def process_segmentation_masks(cell_mask_all_axes,
                               nuclear_mask_all_axes,
                               nucleus_channel,
                               cytoplasm_channel,
                               membrane_channel,
                               image,
                               voxel_size):
	min_JI = 0.0
	max_JI = 0.4
	JI_search_interval = 0.1
	num_steps = int((max_JI - min_JI) / JI_search_interval) + 1
	JI_range = np.linspace(min_JI, max_JI, num_steps)
	print("Matching 2D cells in adjacent slices for each axis...")
	matched_2D_stack_all_JI = {}
	for JI in JI_range:
		matched_2D_stack_all_JI[JI] = {}
		for axis in ['XY', 'XZ', 'YZ']:
			matched_2D_stack_axis = matching_cells_2D(cell_mask_all_axes[axis], JI)
			matched_2D_stack_all_JI[JI][axis] = matched_2D_stack_axis
	
	print("Matching and repairing 3D cells...")
	matched_3D_all_JI = {}
	for JI in JI_range:
		matched_2D_stack_XY = matched_2D_stack_all_JI[JI]['XY']
		matched_2D_stack_XZ = matched_2D_stack_all_JI[JI]['XZ']
		matched_2D_stack_YZ = matched_2D_stack_all_JI[JI]['YZ']
		matched_3D_cell_mask = matching_cells_3D(matched_2D_stack_XY, matched_2D_stack_XZ, matched_2D_stack_YZ)
		matched_3D_all_JI[JI] = matched_3D_cell_mask
	
	print("Matching 3D nuclei...")
	final_matched_3D_cell_mask_JI = {}
	final_matched_3D_nuclear_mask_JI = {}
	for JI in JI_range:
		matched_3D_cell_mask = matched_3D_all_JI[JI]
		final_matched_3D_cell_mask, final_matched_3D_nuclear_mask = matching_nuclei_3D(matched_3D_cell_mask,
		                                                                               nuclear_mask_all_axes['XY'])
		final_matched_3D_cell_mask_JI[JI] = final_matched_3D_cell_mask
		final_matched_3D_nuclear_mask_JI[JI] = final_matched_3D_nuclear_mask
	
	print("Evaluating 3D cell segmentations...")
	quality_score_JI = list()
	metrics_JI = list()
	for JI in JI_range:
		print(f'For JI threshold = {JI}')
		final_matched_3D_cell_mask = final_matched_3D_cell_mask_JI[JI]
		final_matched_3D_nuclear_mask = final_matched_3D_nuclear_mask_JI[JI]
		quality_score, metrics = seg_evaluation_3D(final_matched_3D_cell_mask,
				                                   final_matched_3D_nuclear_mask,
				                                   nucleus_channel,
				                                   cytoplasm_channel,
				                                   membrane_channel,
				                                   image,
				                                   voxel_size,
				                                   './evaluation/model')
		print(f'Quality score = {quality_score}')
		quality_score_JI.append(quality_score)
		metrics_JI.append(metrics)
	
	best_quality_score = max(quality_score_JI)
	best_JI_index = quality_score_JI.index(best_quality_score)
	best_JI = JI_range[best_JI_index]
	best_metrics = metrics_JI[best_JI_index]
	best_cell_mask = final_matched_3D_cell_mask_JI[best_JI]
	best_nuclear_mask = final_matched_3D_nuclear_mask_JI[best_JI]
	
	return best_quality_score, best_metrics, best_cell_mask, best_nuclear_mask


def ThreeDCellComposer(image_path, nucleus_channel_marker_list, cytoplasm_channel_marker_list, membrane_channel_marker_list, segmentation_method, downsample_vector):
	# Process the image
	print("ThreeDCellComposer v1.1")
	print("Generating combined channel files for segmentation...")
	nucleus_channel, cytoplasm_channel, membrane_channel, image = write_IMC_input_channels(image_path,nucleus_channel_marker_list,cytoplasm_channel_marker_list,membrane_channel_marker_list,downsample_vector)
	voxel_size = extract_voxel_size_from_tiff(image_path)
	print('Voxel sizes:',voxel_size)
	
	print("Segmenting every 2D slice across three axes...")
	if segmentation_method in ["deepcell", "custom"]:
		# For a single method
		if segmentation_method == "deepcell":
			cell_mask_all_axes = {}
			nuclear_mask_all_axes = {}
			for axis in ['XY', 'XZ', 'YZ']:
				cell_mask_axis, nuclear_mask_axis = deepcell_segmentation_2D(nucleus_channel, membrane_channel, axis,
				                                                             voxel_size)
				cell_mask_all_axes[axis] = cell_mask_axis
				nuclear_mask_all_axes[axis] = nuclear_mask_axis
		
		
#		elif segmentation_method == "custom":
#			cell_mask_all_axes = {}
#			nuclear_mask_all_axes = {}
#			for axis in ['XY', 'XZ', 'YZ']:
#				cell_mask_axis, nuclear_mask_axis = custom_segmentation(nucleus_channel, cytoplasm_channel,
#				                                                #        membrane_channel, axis, voxel_size)
#				cell_mask_all_axes[axis] = cell_mask_axis
#				nuclear_mask_all_axes[axis] = nuclear_mask_axis
#		
		best_quality_score, best_metrics, best_cell_mask_final, best_nuclear_mask_final = process_segmentation_masks(
			cell_mask_all_axes,
			nuclear_mask_all_axes,
			nucleus_channel,
			cytoplasm_channel,
			membrane_channel,
			image,
			voxel_size)
		
		print(f"Quality Score of final 3D Cell Segmentation = {best_quality_score}")
	
	
#	elif segmentation_method == "compare":
#		# For comparing multiple methods
#		print('installing all methods, it may take some time...')
#		all_methods = ['DeepCell-0.12.6_membrane',
#		               'DeepCell-0.12.6_cytoplasm',
#		               'Cellpose-2.2.2',
#		               'CellProfiler',
#		               'ACSS_classic',
#		               'CellX',
#		               'CellSegm']
#		install_segmentation_methods()
#		split_slices(os.path.dirname(image_path))
#		quality_score_list = list()
#		metrics_list = list()
#		cell_mask_final_list = list()
#		nuclear_mask_final_list = list()
#		for method in all_methods:
#			cell_mask_all_axes, nuclear_mask_all_axes = segmentation_single_method(method, os.path.dirname(image_path), voxel_size)
#			method_quality_score, method_metrics, method_cell_mask_final, method_nuclear_mask_final = process_segmentation_masks(
#				cell_mask_all_axes,
#				nuclear_mask_all_axes,
#				nucleus_channel,
#				cytoplasm_channel,
#				membrane_channel,
#				image,
#				voxel_size)
#			quality_score_list.append(method_quality_score)
#			metrics_list.append(method_metrics)
#			cell_mask_final_list.append(method_cell_mask_final)
#			nuclear_mask_final_list.append(method_nuclear_mask_final)
#
#		best_quality_score = max(quality_score_list)
#		best_quality_score_index = quality_score_list.index(best_quality_score)
#		best_metrics = metrics_list[best_quality_score_index]
#		best_method = all_methods[best_quality_score_index]
#		best_cell_mask_final = cell_mask_final_list[best_quality_score_index]
#		best_nuclear_mask_final = nuclear_mask_final_list[best_quality_score_index]
#		print(f'{best_method} yields the best segmentation.')
#		print(f"Quality Score of final 3D Cell Segmentation = {best_quality_score}")
	
	else:
		print('Invalid segmentation method.')
		exit()

	results_path = f'{os.path.dirname(image_path)}/results'
	if not os.path.exists(results_path):
		os.makedirs(results_path)
	
	mshape= best_cell_mask_final.shape
	#print(mshape,best_nuclear_mask_final.shape)

	import tifffile
	tifffile.imwrite(f'{results_path}/3D_cell_mask.tif', best_cell_mask_final)
	tifffile.imwrite(f'{results_path}/3D_nuclear_mask.tif', best_nuclear_mask_final)
	combined_mask = np.zeros([2,mshape[0],mshape[1],mshape[2]], dtype='uint16')
	combined_mask[0,:,:,:]=best_cell_mask_final
	combined_mask[1,:,:,:]=best_nuclear_mask_final
	#print(combined_mask.shape)

	swapped = np.swapaxes(combined_mask,0,1)
	#print(swapped.shape)
	#swapped = np.swapaxes(swapped,2,3)
	#print(swapped.shape)

	image_labels = [f'{i}' for i in range(combined_mask.shape[0] * combined_mask.shape[1])]
	tifffile.imwrite(f'{results_path}/3D_cell_nuclear_mask.tif', swapped,
		imagej=True,
		         resolution=(1./float(voxel_size[0]), 1./float(voxel_size[1])),
		metadata={
			'spacing': float(voxel_size[2]),
			'unit': 'um',
			'finterval': 1/10,
			'fps': 10.0,
			'axes': 'ZCYX',
			'Labels': image_labels,
			'Image': {
				'Pixels': {
					'PhysicalSizeX': voxel_size[0],
					#'PhysicalSizeXUnit': "µm",
 					'PhysicalSizeY': voxel_size[1], 
 					#'PhysicalSizeYUnit': "µm",
 					'PhysicalSizeZ': voxel_size[2], 
					#'PhysicalSizeZUnit': "µm",
 				},
 				'Channel': { 'ID': 1, 'Name': "Cell",},
				'Channel': { 'ID': 2, 'Name': "Nucleus",},
			},
		}
	)

	metrics_path = f'{results_path}/metrics.json'
	
	with open(metrics_path, 'w') as f:
		json.dump(best_metrics, f)
	
	quality_score_path = f'{results_path}/quality_score.txt'
	np.savetxt(quality_score_path, [best_quality_score], fmt='%f')

	print("Generating surface meshes for visualization in Blender..")
	best_cell_mask_final_colored, number_of_colors = coloring_3D(best_cell_mask_final)
	meshing_3D(best_cell_mask_final, best_cell_mask_final_colored, number_of_colors, results_path)
	
	print("3D Segmentation and Evaluation Completed.")
