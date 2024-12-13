import argparse
from pathlib import Path
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
import json
import os
import pickle
import numpy as np
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
from scipy.ndimage import binary_erosion, binary_dilation
import tifffile

from ome_utils import find_ome_tiffs
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

"""
3DCellComposer MAIN PROGRAM
Author: Haoran Chen and Robert F. Murphy
Version: 1.1 December 14, 2023 R.F.Murphy, Haoran Chen
        Fixes in
        deepcell_only.py
        single_method_eval_3D.py
        meshing_3D.py
"""

def parse_marker_list(arg):
    return arg.split(',')

def convert_to_sequential_labels_parallel(mask):
    """Parallelized version of sequential label conversion"""
    unique_labels = np.unique(mask)
    unique_labels = unique_labels[unique_labels > 0]
    
    # Create relabeling map
    relabel_map = {old_label: new_label for new_label, old_label 
                   in enumerate(unique_labels, start=1)}
    
    # Process chunks in parallel
    def process_chunk(chunk_slice):
        chunk = mask[chunk_slice]
        result = np.zeros_like(chunk, dtype=np.int32)
        for old_label, new_label in relabel_map.items():
            result[chunk == old_label] = new_label
        return chunk_slice, result

    # Split into chunks based on CPU count
    num_cores = mp.cpu_count()
    chunks = np.array_split(range(mask.shape[0]), num_cores)
    chunk_slices = [slice(chunk[0], chunk[-1] + 1) for chunk in chunks]
    
    sequential_mask = np.zeros_like(mask, dtype=np.int32)
    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_chunk, chunk_slice) 
                  for chunk_slice in chunk_slices]
        
        for future in futures:
            chunk_slice, result = future.result()
            sequential_mask[chunk_slice] = result
    
    return sequential_mask, relabel_map

def process_slice_boundaries(args):
    """Process a single slice for boundary detection"""
    z_idx, mask_slice = args
    boundaries = np.zeros_like(mask_slice)
    labels = np.unique(mask_slice)[1:]  # exclude 0
    
    for label in labels:
        binary_mask = (mask_slice == label)
        eroded = binary_erosion(binary_mask)
        boundary = binary_mask & ~eroded
        boundaries[boundary] = label
    
    return z_idx, boundaries

def get_3D_boundaries_parallel(mask):
    """Parallelized version of 3D boundary detection"""
    sequential_mask, relabel_map = convert_to_sequential_labels_parallel(mask)
    print(f"Relabeled {len(relabel_map)} objects")
    
    boundaries = np.zeros_like(sequential_mask)
    
    # Prepare slices for parallel processing
    slice_args = [(z, sequential_mask[z]) for z in range(sequential_mask.shape[0])]
    
    # Process in parallel using ProcessPoolExecutor
    num_cores = mp.cpu_count()
    with ProcessPoolExecutor(max_workers=num_cores) as executor:
        futures = list(executor.map(process_slice_boundaries, slice_args))
        
        for z_idx, slice_boundaries in futures:
            boundaries[z_idx] = slice_boundaries
            if z_idx % 10 == 0:
                print(f"Processed slice {z_idx}/{sequential_mask.shape[0]}")
    
    return sequential_mask, boundaries

def parallel_matching_2D(cell_mask_all_axes, JI):
    """Parallelized version of 2D cell matching"""
    matched_2D_stack_all_axes = {}
    
    with ProcessPoolExecutor() as executor:
        futures = {
            axis: executor.submit(matching_cells_2D, cell_mask_all_axes[axis], JI)
            for axis in ['XY', 'XZ', 'YZ']
        }
        
        for axis, future in futures.items():
            matched_2D_stack_all_axes[axis] = future.result()
            
    return matched_2D_stack_all_axes

def process_segmentation_masks_parallel(cell_mask_all_axes,
                                      nuclear_mask_all_axes,
                                      nucleus_channel,
                                      cytoplasm_channel,
                                      membrane_channel,
                                      image,
                                      voxel_size):
    """Parallelized version of segmentation mask processing"""
    min_JI = 0.0
    max_JI = 0.4
    JI_search_interval = 0.1
    JI_range = np.linspace(min_JI, max_JI, int((max_JI - min_JI) / JI_search_interval) + 1)
    
    print("Matching 2D cells in parallel across all axes...")
    with ProcessPoolExecutor() as executor:
        futures = {
            JI: executor.submit(parallel_matching_2D, cell_mask_all_axes, JI)
            for JI in JI_range
        }
        matched_2D_stack_all_JI = {JI: future.result() 
                                  for JI, future in futures.items()}
    
    print("Matching and repairing 3D cells in parallel...")
    with ProcessPoolExecutor() as executor:
        futures = {
            JI: executor.submit(
                matching_cells_3D,
                matched_2D_stack_all_JI[JI]['XY'],
                matched_2D_stack_all_JI[JI]['XZ'],
                matched_2D_stack_all_JI[JI]['YZ']
            )
            for JI in JI_range
        }
        matched_3D_all_JI = {JI: future.result() 
                            for JI, future in futures.items()}
    
    print("Matching 3D nuclei and evaluating in parallel...")
    def process_JI(JI):
        matched_3D_cell_mask = matched_3D_all_JI[JI]
        final_matched_3D_cell_mask, final_matched_3D_nuclear_mask = matching_nuclei_3D(
            matched_3D_cell_mask, nuclear_mask_all_axes['XY']
        )
        
        quality_score, metrics = seg_evaluation_3D(
            final_matched_3D_cell_mask,
            final_matched_3D_nuclear_mask,
            nucleus_channel,
            cytoplasm_channel,
            membrane_channel,
            image,
            voxel_size,
            './evaluation/model'
        )
        
        return JI, quality_score, metrics, final_matched_3D_cell_mask, final_matched_3D_nuclear_mask
    
    with ProcessPoolExecutor() as executor:
        results = list(executor.map(process_JI, JI_range))
    
    # Find best result
    best_result = max(results, key=lambda x: x[1])  # x[1] is quality_score
    _, best_quality_score, best_metrics, best_cell_mask, best_nuclear_mask = best_result
    
    return best_quality_score, best_metrics, best_cell_mask, best_nuclear_mask

def parallel_segmentation(nucleus_channel, membrane_channel, voxel_size):
    """Parallelized version of axis segmentation"""
    cell_mask_all_axes = {}
    nuclear_mask_all_axes = {}
    
    with ProcessPoolExecutor() as executor:
        # Process all axes in parallel
        futures = [
            executor.submit(deepcell_segmentation_2D, nucleus_channel, membrane_channel, axis, voxel_size)
            for axis in ['XY', 'XZ', 'YZ']
        ]
        
        for axis, future in zip(['XY', 'XZ', 'YZ'], futures):
            cell_mask_axis, nuclear_mask_axis = future.result()
            cell_mask_all_axes[axis] = cell_mask_axis
            nuclear_mask_all_axes[axis] = nuclear_mask_axis
            
    return cell_mask_all_axes, nuclear_mask_all_axes

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
    print("3DCellComposer v1.1")
    print("Generating input channels for segmentation...")
    nucleus_channel, cytoplasm_channel, membrane_channel, image = write_IMC_input_channels(
        image_path,
        args.results_path,
        args.nucleus_channel_marker_list,
        args.cytoplasm_channel_marker_list,
        args.membrane_channel_marker_list
    )
    voxel_size = extract_voxel_size_from_tiff(image_path)
    
    print("Segmenting every 2D slices across three axes...")
    if args.segmentation_method in ["deepcell", "custom"]:
        if args.segmentation_method == "deepcell":
            # Check for cached results
            cell_mask_file_name = 'cell_mask_all_axes.pkl'
            nuclear_mask_file_name = 'nuclear_mask_all_axes.pkl'

            if os.path.exists(cell_mask_file_name) and os.path.exists(nuclear_mask_file_name):
                print("Loading cached segmentation results...")
                with open(cell_mask_file_name, 'rb') as f:
                    cell_mask_all_axes = pickle.load(f)
                with open(nuclear_mask_file_name, 'rb') as f:
                    nuclear_mask_all_axes = pickle.load(f)
            else:
                print("Running parallel segmentation...")
                cell_mask_all_axes, nuclear_mask_all_axes = parallel_segmentation(
                    nucleus_channel, membrane_channel, voxel_size
                )
                
                # Cache results
                with open(cell_mask_file_name, 'wb') as f:
                    pickle.dump(cell_mask_all_axes, f)
                with open(nuclear_mask_file_name, 'wb') as f:
                    pickle.dump(nuclear_mask_all_axes, f)
                    
        elif args.segmentation_method == "custom":
            cell_mask_all_axes, nuclear_mask_all_axes = parallel_segmentation(
                nucleus_channel, cytoplasm_channel, membrane_channel, voxel_size,
                segmentation_func=custom_segmentation
            )

        # Process masks in parallel
        best_quality_score, best_metrics, best_cell_mask_final, best_nuclear_mask_final = \
            process_segmentation_masks_parallel(
                cell_mask_all_axes,
                nuclear_mask_all_axes,
                nucleus_channel,
                cytoplasm_channel,
                membrane_channel,
                image,
                voxel_size
            )
        
        print(f"Quality Score of this 3D Cell Segmentation = {best_quality_score}")
    
    elif args.segmentation_method == "compare":
        print('Installing all methods, it may take some time...')
        all_methods = [
            'DeepCell-0.12.6_membrane',
            'DeepCell-0.12.6_cytoplasm',
            'Cellpose-2.2.2',
            'CellProfiler',
            'ACSS_classic',
            'CellX',
            'CellSegm'
        ]
        install_segmentation_methods()
        split_slices(image_path.parent)
        
        # Process all methods in parallel
        def process_method(method):
            cell_mask_all_axes, nuclear_mask_all_axes = segmentation_single_method(
                method, image_path.parent, voxel_size
            )
            return process_segmentation_masks_parallel(
                cell_mask_all_axes,
                nuclear_mask_all_axes,
                nucleus_channel,
                cytoplasm_channel,
                membrane_channel,
                image,
                voxel_size
            )
        
        with ProcessPoolExecutor() as executor:
            method_results = list(executor.map(process_method, all_methods))
        
        # Find best method
        quality_scores = [result[0] for result in method_results]
        best_idx = np.argmax(quality_scores)
        best_quality_score, best_metrics, best_cell_mask_final, best_nuclear_mask_final = method_results[best_idx]
        best_method = all_methods[best_idx]
        
        print(f'{best_method} yields the best segmentation.')
        print(f"Quality Score of this 3D Cell Segmentation = {best_quality_score}")

    # Generate and save boundaries in parallel
    print("Generating boundaries...")
    cell_mask_sequential, cell_boundaries = get_3D_boundaries_parallel(best_cell_mask_final)
    nuclear_mask_sequential, nuclear_boundaries = get_3D_boundaries_parallel(best_nuclear_mask_final)

    # Write results
    print("Saving results...")
    tifffile.imwrite(args.results_path / '3D_cell_mask.tif', cell_mask_sequential)
    tifffile.imwrite(args.results_path / '3D_nuclear_mask.tif', nuclear_mask_sequential)
    tifffile.imwrite(args.results_path / '3D_cell_boundaries.tif', cell_boundaries)
    tifffile.imwrite(args.results_path / '3D_nuclear_boundaries.tif', nuclear_boundaries)

    with open(args.results_path / 'metrics.json', 'w') as f:
        json.dump(best_metrics, f)

    np.savetxt(args.results_path / 'quality_score.txt', [best_quality_score], fmt='%f')

    print("Generating surface meshes for visualization in Blender...")
    best_cell_mask_final_colored, number_of_colors = coloring_3D(best_cell_mask_final)
    meshing_3D(best_cell_mask_final, best_cell_mask_final_colored, number_of_colors, args.results_path)
    
    print("3D Segmentation and Evaluation Completed.")

if __name__ == "__main__":
    main()