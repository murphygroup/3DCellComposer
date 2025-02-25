import argparse
from pathlib import Path
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
import json
from datetime import datetime

from ome_utils import find_ome_tiffs
import pickle

import tifffile

from segmentation_2D.preprocessing import *
from segmentation_2D.deepcell_only import deepcell_segmentation_2D
from segmentation_2D.wrapper.custom_segmentation_wrapper import custom_segmentation
from segmentation_2D.installation.install_all_methods import install_segmentation_methods
from segmentation_2D.all_methods_segmentation import *
from segmentation_3D.match_2D_cells import matching_cells_2D
from segmentation_3D.match_3D_cells import matching_cells_3D
from segmentation_3D.match_3D_nuclei import matching_nuclei_3D
from evaluation.single_method_eval_3D import seg_evaluation_3D
from visualization.coloring_3D import coloring_3D
from visualization.meshing_3D import meshing_3D

from scipy.ndimage import binary_dilation, binary_erosion

"""
3DCellComposer MAIN PROGRAM
Author: Haoran Chen and Robert F. Murphy
Version: 1.1 December 14, 2023 R.F.Murphy, Haoran Chen
        Fixes in
        deepcell_only.py
        single_method_eval_3D.py
        meshing_3D.py
Version: 1.3 February 15, 2025 R.F.Murphy, Ted Zhang
        Correct sampling code to use replication to fill intermediate slices
        Allow different sampling intervals for XY, XZ and YZ
        Convert to even spacing by dividing by specific factor
        Save and reuse previous segmented slices
        Optionally save results after each trial
"""

def parse_marker_list(arg):
	return arg.split(',')


def process_segmentation_masks(cell_mask_all_axes,
                               nuclear_mask_all_axes,
                               nucleus_channel,
                               cytoplasm_channel,
                               membrane_channel,
                               image,
                               voxel_size,max_JI):
	min_JI = 0.0
	#max_JI = 0.4
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
			matched_2D_stack_all_JI[JI][axis] = matched_2D_stack_axis
			matched_2D_stack_all_JI[JI][axis] = matched_2D_stack_axis
	print(f"{datetime.now()} Matching and repairing 3D cells...")
	matched_3D_all_JI = {}
	for JI in JI_range:
		matched_2D_stack_XY = matched_2D_stack_all_JI[JI]['XY']
		matched_2D_stack_XZ = matched_2D_stack_all_JI[JI]['XZ']
		matched_2D_stack_YZ = matched_2D_stack_all_JI[JI]['YZ']
		matched_3D_cell_mask = matching_cells_3D(matched_2D_stack_XY, matched_2D_stack_XZ, matched_2D_stack_YZ)
		matched_3D_all_JI[JI] = matched_3D_cell_mask

	print(f"{datetime.now()} Matching 3D nuclei...")
	final_matched_3D_cell_mask_JI = {}
	final_matched_3D_nuclear_mask_JI = {}
	for JI in JI_range:
		matched_3D_cell_mask = matched_3D_all_JI[JI]
		final_matched_3D_cell_mask, final_matched_3D_nuclear_mask = matching_nuclei_3D(matched_3D_cell_mask,
		                                                                               nuclear_mask_all_axes['XY'])
		#print(final_matched_3D_cell_mask.shape,final_matched_3D_nuclear_mask.shape)
		#final_matched_3D_cell_mask = fill_in_slices3D(final_matched_3D_cell_mask, nucleus_channel.shape)
		#final_matched_3D_nuclear_mask = fill_in_slices3D(final_matched_3D_nuclear_mask, nucleus_channel.shape)
		#print(final_matched_3D_cell_mask.shape,final_matched_3D_nuclear_mask.shape)
		final_matched_3D_cell_mask_JI[JI] = final_matched_3D_cell_mask
		final_matched_3D_nuclear_mask_JI[JI] = final_matched_3D_nuclear_mask
	
	print("{datetime.now()} Evaluating 3D cell segmentation...")
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
		print(f'- quality score = {quality_score}')
		quality_score_JI.append(quality_score)
		metrics_JI.append(metrics)
	
	best_quality_score = max(quality_score_JI)
	best_JI_index = quality_score_JI.index(best_quality_score)
	best_JI = JI_range[best_JI_index]
	best_metrics = metrics_JI[best_JI_index]
	best_cell_mask = final_matched_3D_cell_mask_JI[best_JI]
	best_nuclear_mask = final_matched_3D_nuclear_mask_JI[best_JI]
	print(f"{datetime.now()}")
	
	return best_quality_score, best_metrics, best_cell_mask, best_nuclear_mask

def convert_to_sequential_labels(mask):

    unique_labels = np.unique(mask)
    unique_labels = unique_labels[unique_labels > 0]
    
    # Create relabeling map (old label -> new sequential label)
    relabel_map = {old_label: new_label for new_label, old_label 
                   in enumerate(unique_labels, start=1)}
    
    # Create new sequential mask
    sequential_mask = np.zeros_like(mask, dtype=np.int32)
    for old_label, new_label in relabel_map.items():
        sequential_mask[mask == old_label] = new_label
    
    return sequential_mask, relabel_map

def get_3D_boundaries(mask):

    # First relabel the mask with sequential integers
    sequential_mask, relabel_map = convert_to_sequential_labels(mask)
    print(f"Relabeled {len(relabel_map)} objects")
    
    # Generate boundaries with sequential labels
    boundaries = np.zeros_like(sequential_mask)
    # for label in range(1, len(relabel_map) + 1):
    #     binary_mask = (sequential_mask == label)
    #     dilated = binary_dilation(binary_mask)
    #     boundaries[dilated & ~binary_mask] = label
    # Process each z-slice
    for z in range(sequential_mask.shape[0]):
        boundaries[z] = get_2D_boundaries(sequential_mask[z])
        
        # Print progress every 10 slices
        if z % 10 == 0:
            print(f"Processed slice {z}/{sequential_mask.shape[0]}")
    
    return sequential_mask, boundaries

def get_2D_boundaries(mask_slice):

    boundaries = np.zeros_like(mask_slice)
    labels = np.unique(mask_slice)[1:]  # exclude 0
    
    for label in labels:
        # Create binary mask for current label
        binary_mask = (mask_slice == label)
        
        # Get boundary by subtracting eroded mask from original
        eroded = binary_erosion(binary_mask)
        boundary = binary_mask & ~eroded
        
        # Set boundary with original label value
        boundaries[boundary] = label
    
    return boundaries

def writeresults(rpath,best_cell_mask_final,best_nuclear_mask_final,best_metrics,best_quality_score):
	if not rpath.is_dir():
		rpath.mkdir(exist_ok=True, parents=True)
	#get boundaries
	cell_mask, cell_boundaries = get_3D_boundaries(best_cell_mask_final)
	nuclear_mask, nuclear_boundaries = get_3D_boundaries(best_nuclear_mask_final)
	# Write masks
	tifffile.imwrite(rpath / '3D_cell_mask.tif', cell_mask)
	tifffile.imwrite(rpath / '3D_nuclear_mask.tif', nuclear_mask)
	tifffile.imwrite(rpath / '3D_cell_boundaries.tif', cell_boundaries)
	tifffile.imwrite(rpath / '3D_nuclear_boundaries.tif', nuclear_boundaries)

	with open(rpath / 'metrics.json', 'w') as f:
		json.dump(best_metrics, f)

	np.savetxt(rpath / 'quality_score.txt', [best_quality_score], fmt='%f')

def main():
	parser = argparse.ArgumentParser(description="3DCellComposer Image Processor")
	parser.add_argument("image_path",
						help="Path to the multiplexed image, or directory containing one image",
						type=Path)
	parser.add_argument("nucleus_channel_marker_list", type=parse_marker_list,
	                    help="A list of nuclear marker(s) in multiplexed image as input for segmentation")
	parser.add_argument("cytoplasm_channel_marker_list", type=parse_marker_list,
	                    help="A list of cytoplasmic marker(s) in multiplexed image as input for segmentation")
	parser.add_argument("membrane_channel_marker_list", type=parse_marker_list,
	                    help="A list of cell membrane marker(s) in multiplexed image as input for segmentation")
	parser.add_argument("--segmentation_method", type=str, choices=["deepcell", "compare", "custom"],
	                    default="deepcell",
	                    help="Choose 2D segmentation method: 'deepcell' for DeepCell segmentation which performed the best in our evaluation, 'compare' to compare and select the best among 7 methods, or 'custom' to use a user-provided method.")
	parser.add_argument('--results_path', type=Path, default=Path('results'),
						help="Path where results will be written")
	parser.add_argument('--chunk_size', type=int, default=100,
						help="Chunk size for segmentation")
	parser.add_argument('--sampling_interval', type=parse_marker_list, default='1,1,1',
						help="Sampling interval for segmentation in XY, XZ and YZ directions")
	#parser.add_argument('--interpolation_method', type=str, default='cubic',
	#					help="Interpolation method for segmentation")
	#parser.add_argument('--fill_value', type=str, default='extrapolate',
	#					help="Fill value for segmentation")
	parser.add_argument('--max_tries', type=int, default=10,
						help="Maximum number of tries for sampled segmentation")
	parser.add_argument('--quality_threshold', type=float, default=float("inf"),
						help="Quality threshold for segmentation")
	parser.add_argument('--sampling_reduce', type=int, default=2,
						help="Divisor to reduce sampling interval for segmentation")
	parser.add_argument('--max_JI', type=float, default=0.4,
						help="Upper limit for Jaccard index for cell merging")




	args = parser.parse_args()

	if args.image_path.is_dir():
		ome_tiffs = list(find_ome_tiffs(args.image_path))
		image_path = ome_tiffs[0]
	elif args.image_path.is_file():
		image_path = args.image_path
	else:
		raise FileNotFoundError(f'Not a file or directory: {args.image_path}')

	if not args.results_path.is_dir():
		args.results_path.mkdir(exist_ok=True, parents=True)

	# Process the image
	print("3DCellComposer v1.3")
	print("Generating input channels for segmentation...")
	nucleus_channel, cytoplasm_channel, membrane_channel, image = write_IMC_input_channels(image_path,
																						   args.results_path,
	                                                                                       args.nucleus_channel_marker_list,
	                                                                                       args.cytoplasm_channel_marker_list,
	                                                                                       args.membrane_channel_marker_list)
	#print(nucleus_channel.shape,cytoplasm_channel.shape,membrane_channel.shape,image.shape)
	print(f"Marker channels shape: {nucleus_channel.shape}, All channels shape: {image.shape}")
	voxel_size = extract_voxel_size_from_tiff(image_path)
	print(voxel_size)

	#segmentation
	print(type(args.sampling_interval))
	print(args.sampling_interval)
	svals = args.sampling_interval
	sampling_interval = {'XY': int(svals[0]), 'XZ': int(svals[1]), 'YZ': int(svals[2])}
	print(f"Sampling interval {sampling_interval}")
	max_tries = args.max_tries
	quality_threshold = args.quality_threshold
	sampling_reduce = args.sampling_reduce

	print("Segmenting every 2D slices across three axes...")
	if args.segmentation_method in ["deepcell", "custom"]:

		while max_tries > 0:
			# For a single method

			if args.segmentation_method == "deepcell":


				cell_mask_all_axes = {}
				nuclear_mask_all_axes = {}
				for axis in ['XY', 'XZ', 'YZ']:
					cell_mask_axis, nuclear_mask_axis = deepcell_segmentation_2D(nucleus_channel, membrane_channel, axis,
																				voxel_size, sampling_interval[axis], args.chunk_size, args.results_path)
					cell_mask_all_axes[axis] = cell_mask_axis
					nuclear_mask_all_axes[axis] = nuclear_mask_axis

					#save for future debug
					#save dict
					#with open('cell_mask_all_axes.pkl', 'wb') as cell_mask_file:
					#	pickle.dump(cell_mask_all_axes, cell_mask_file)

					#with open('nuclear_mask_all_axes.pkl', 'wb') as nuclear_mask_file:
					#	pickle.dump(nuclear_mask_all_axes, nuclear_mask_file)

			elif args.segmentation_method == "custom":
				cell_mask_all_axes = {}
				nuclear_mask_all_axes = {}
				for axis in ['XY', 'XZ', 'YZ']:
					cell_mask_axis, nuclear_mask_axis = custom_segmentation(nucleus_channel, cytoplasm_channel,
																			membrane_channel, axis, voxel_size)
					cell_mask_all_axes[axis] = cell_mask_axis
					nuclear_mask_all_axes[axis] = nuclear_mask_axis
			best_quality_score, best_metrics, best_cell_mask_final, best_nuclear_mask_final = process_segmentation_masks(
				cell_mask_all_axes,
				nuclear_mask_all_axes,
				nucleus_channel,
				cytoplasm_channel,
				membrane_channel,
				image,
				voxel_size, args.max_JI)
		
			print(f"Quality Score of this 3D Cell Segmentation = {best_quality_score}")
			trialpath = args.results_path / ('trial' + str(args.max_tries-max_tries))
			print(trialpath)
			writeresults(trialpath,best_cell_mask_final,best_nuclear_mask_final,best_metrics,best_quality_score)
			if best_quality_score > args.quality_threshold or max(sampling_interval.values())==1 or max_tries ==1:
				break
			else:
				max_tries -= 1
				for axis in ['XY', 'XZ', 'YZ']:
					sampling_interval[axis] = int(max(sampling_interval[axis] / sampling_reduce,1))
				print(f"Quality score is too low, Sampling interval is reduced to {sampling_interval}")
				print(f"Tries left: {max_tries}")

	

	
	elif args.segmentation_method == "compare":

		# For comparing multiple methods
		print('installing all methods, it may take some time...')
		all_methods = ['DeepCell-0.12.6_membrane',
		               'DeepCell-0.12.6_cytoplasm',
		               'Cellpose-2.2.2',
		               'CellProfiler',
		               'ACSS_classic',
		               'CellX',
		               'CellSegm']
		install_segmentation_methods()
		split_slices(image_path.parent)
		quality_score_list = list()
		metrics_list = list()
		cell_mask_final_list = list()
		nuclear_mask_final_list = list()
		for method in all_methods:
			cell_mask_all_axes, nuclear_mask_all_axes = segmentation_single_method(method, image_path.parent, voxel_size)
			method_quality_score, method_metrics, method_cell_mask_final, method_nuclear_mask_final = process_segmentation_masks(
				cell_mask_all_axes,
				nuclear_mask_all_axes,
				nucleus_channel,
				cytoplasm_channel,
				membrane_channel,
				image,
				voxel_size)
			quality_score_list.append(method_quality_score)
			metrics_list.append(method_metrics)
			cell_mask_final_list.append(method_cell_mask_final)
			nuclear_mask_final_list.append(method_nuclear_mask_final)
		
		best_quality_score = max(quality_score_list)
		best_quality_score_index = quality_score_list.index(best_quality_score)
		best_metrics = metrics_list[best_quality_score_index]
		best_method = all_methods[best_quality_score_index]
		best_cell_mask_final = cell_mask_final_list[best_quality_score_index]
		best_nuclear_mask_final = nuclear_mask_final_list[best_quality_score_index]
		print(f'{best_method} yields the best segmentation.')
		print(f"Quality Score of this 3D Cell Segmentation = {best_quality_score}")

	writeresults(args.results_path,best_cell_mask_final,best_nuclear_mask_final,best_metrics,best_quality_score)

	print("Generating surface meshes for visualization in Blender..")
	best_cell_mask_final_colored, number_of_colors = coloring_3D(best_cell_mask_final)
	meshing_3D(best_cell_mask_final, best_cell_mask_final_colored, number_of_colors, args.results_path)
	
	print("3D Segmentation and Evaluation Completed.")


if __name__ == "__main__":
	main()
