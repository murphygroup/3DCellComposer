import argparse
import os.path
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
import json
from ThreeDCellComposer.ThreeDCellComposer import ThreeDCellComposer

"""
3DCellComposer MAIN PROGRAM
Author: Haoran Chen and Robert F. Murphy
Version: 1.1 December 14, 2023 R.F.Murphy, Haoran Chen
        Fixes in
        deepcell_only.py
        single_method_eval_3D.py
        meshing_3D.py
         1.2 February 7, 2024 R.F.Murphy
        create and use PyPi package
Version: 1.3 February 15, 2025 R.F.Murphy, Ted Zhang
        Correct sampling code to use replication to fill intermediate slices
        Allow different sampling intervals for XY, XZ and YZ
        Convert to even spacing by dividing by specific factor
        Save and reuse previous segmented slices
        Optionally save results after each trial
Version: 1.4 February 20, 2025 R.F.Murphy
        Allow downsampling before cell segmentation and matching (and upsampling after)
        Allow skipping of evaluation and/or creating blender output
Version: 1.5 March 7, 2025 R.F.Murphy
        Save command line arguments to results folder
        Estimate completion time
"""

def parse_marker_list(arg):
	return arg.split(',')

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
	parser.add_argument('--JI_range', type=parse_marker_list, default='0.0,0.4,5',
						help="Range for Jaccard index for cell merging")
	parser.add_argument('--skip_eval', type=bool, default=False,
						help="Skip CellSegmentationEvaluator")
	parser.add_argument('--skip_blender', type=bool, default=False,
						help="Skip generating blender files")
	parser.add_argument('--downsample_vector', type=parse_marker_list, default="1,1,1",
						help="Vector for downsampling each axis before segmentation")
	parser.add_argument('--maxima_threshold', type=float, default=0.075,
						help="DeepCell parameter for segmentation (lower gives more cells)")
	parser.add_argument('--interior_threshold', type=float, default=0.2,
						help="DeepCell parameter for segmentation (lower gives fewer cells)")
	parser.add_argument('--compartment', type=str, default="both",
						help="DeepCell channels to use (both, whole-cell, nuclear)")
	parser.add_argument('--crop_limits', type=parse_marker_list, default="0,-1,0,-1,0,-1",
						help="Zl,Zh,Yl,Yh,Xl,Xh limits for cropping before segmentation")
	parser.add_argument('--min_slice_padding', type=int, default="512",
						help="minimum size to pad slices to")
	
	args = parser.parse_args()
	ThreeDCellComposer(args.image_path,args.nucleus_channel_marker_list,args.cytoplasm_channel_marker_list,args.membrane_channel_marker_list,args.segmentation_method,args.results_path,args.chunk_size,args.sampling_interval,args.max_tries,args.quality_threshold,args.sampling_reduce,args.JI_range,args.skip_eval,args.skip_blender,args.downsample_vector,args.maxima_threshold,args.interior_threshold,args.compartment,args.crop_limits,args.min_slice_padding)
	print("3DCellComposer completed.")


if __name__ == "__main__":
	main()
