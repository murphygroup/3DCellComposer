import argparse
import os.path
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
import json

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
	max_JI = 0.5
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
	
	print("Evaluating 3D cell segmentation...")
	quality_score_JI = list()
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
	
	return best_quality_score, best_metrics, best_cell_mask, best_nuclear_mask


def main():
	parser = argparse.ArgumentParser(description="3DCellComposer Image Processor")
	parser.add_argument("image_path", help="Path to the multiplexed image")
	parser.add_argument("nucleus_channel_marker_list", type=parse_marker_list,
	                    help="A list of nuclear marker(s) in multiplexed image as input for segmentation")
	parser.add_argument("cytoplasm_channel_marker_list", type=parse_marker_list,
	                    help="A list of cytoplasmic marker(s) in multiplexed image as input for segmentation")
	parser.add_argument("membrane_channel_marker_list", type=parse_marker_list,
	                    help="A list of cell membrane marker(s) in multiplexed image as input for segmentation")
	parser.add_argument("--segmentation_method", type=str, choices=["deepcell", "compare", "custom"],
	                    default="deepcell",
	                    help="Choose 2D segmentation method: 'deepcell' for DeepCell segmentation which performed the best in our evaluation, 'compare' to compare and select the best among 7 methods, or 'custom' to use a user-provided method.")
	
	args = parser.parse_args()
	# Process the image
	print("Generating input channels for segmentation...")
	nucleus_channel, cytoplasm_channel, membrane_channel, image = write_IMC_input_channels(args.image_path,
	                                                                                       args.nucleus_channel_marker_list,
	                                                                                       args.cytoplasm_channel_marker_list,
	                                                                                       args.membrane_channel_marker_list)
	voxel_size = extract_voxel_size_from_tiff(args.image_path)
	
	print("Segmenting every 2D slices across three axes...")
	if args.segmentation_method in ["deepcell", "custom"]:
		# For a single method
		if args.segmentation_method == "deepcell":
			cell_mask_all_axes = {}
			nuclear_mask_all_axes = {}
			for axis in ['XY', 'XZ', 'YZ']:
				cell_mask_axis, nuclear_mask_axis = deepcell_segmentation_2D(nucleus_channel, membrane_channel, axis,
				                                                             voxel_size)
				cell_mask_all_axes[axis] = cell_mask_axis
				nuclear_mask_all_axes[axis] = nuclear_mask_axis
		
		
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
			voxel_size)
		
		print(f"Quality Score of this 3D Cell Segmentation = {quality_score}")
	
	
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
		split_slices(os.path.dirname(args.image_path))
		quality_score_list = list()
		metrics_list = list()
		cell_mask_final_list = list()
		nuclear_mask_final_list = list()
		for method in all_methods:
			cell_mask_all_axes, nuclear_mask_all_axes = segmentation_single_method(method, os.path.dirname(
				"./data/3D_IMC_image.ome.tiff"), voxel_size)
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
	
	results_path = f'{os.path.dirname(args.image_path)}/results'
	if not os.path.exists(results_path):
		os.makedirs(results_path)
	
	import tifffile
	tifffile.imwrite(f'{results_path}/3D_cell_mask.tif', best_cell_mask_final)
	tifffile.imwrite(f'{results_path}/3D_nuclear_mask.tif', best_nuclear_mask_final)
	
	metrics_path = f'{results_path}/metrics.json'
	
	with open(metrics_path, 'w') as f:
		json.dump(best_metrics, f)
	
	quality_score_path = f'{results_path}/quality_score.txt'
	np.savetxt(quality_score_path, [best_quality_score], fmt='%f')

	print("Generating surface meshes for visualization in Blender..")
	best_cell_mask_final_colored, number_of_colors = coloring_3D(best_cell_mask_final)
	meshing_3D(best_cell_mask_final, best_cell_mask_final_colored, number_of_colors)
	
	print("3D Segmentation and Evaluation Completed.")


if __name__ == "__main__":
	main()
