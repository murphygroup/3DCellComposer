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
	
	args = parser.parse_args()
	ThreeDCellComposer(args.image_path,args.nucleus_channel_marker_list,args.cytoplasm_channel_marker_list,args.membrane_channel_marker_list,args.segmentation_method)
	print("3DCellComposer completed.")


if __name__ == "__main__":
	main()
