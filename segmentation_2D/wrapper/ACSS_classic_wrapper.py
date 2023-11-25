import os
import sys
from skimage.io import imread
from os.path import join
from skimage.measure import label
import numpy as np
import matplotlib.pyplot as plt
from aicssegmentation.core.seg_dot import dot_2d_slice_by_slice_wrapper
from aicssegmentation.core.pre_processing_utils import intensity_normalization, image_smoothing_gaussian_3d
from skimage.morphology import remove_small_objects, binary_closing, ball, disk, erosion, dilation   # function for post-processing (size filter)
from aicssegmentation.core.MO_threshold import MO
from skimage.io import imsave
import bz2
import pickle

os.getcwd()

file_dir = sys.argv[1]
Data1 = imread(join(file_dir, 'membrane.tif'))
Data2 = imread(join(file_dir, 'nucleus.tif'))

downsample_percentage = int(sys.argv[2])
pixel_size_in_nano = int(sys.argv[3])

Data1 = np.expand_dims(Data1, axis=0)
Data2 = np.expand_dims(Data2, axis=0)


struct_img0 = Data1.copy()
struct_img1 = Data2.copy()


intensity_scaling_param = [0.5, 15]
gaussian_smoothing_sigma = 1

struct_img0 = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)
struct_img1 = intensity_normalization(struct_img1, scaling_param=intensity_scaling_param)

struct_img = struct_img0 + struct_img1
structure_img_smooth = image_smoothing_gaussian_3d(struct_img, sigma=gaussian_smoothing_sigma)

minArea = 5 * (1000 ** 2) / (pixel_size_in_nano ** 2) * (downsample_percentage ** 2) / (100 ** 2)


bw, object_for_debug = MO(structure_img_smooth, global_thresh_method='ave', object_minArea=minArea, return_object=True)
s2_param_bright = [[2, 0.025]]
s2_param_dark = [[2, 0.025], [1, 0.025]]

bw_extra = dot_2d_slice_by_slice_wrapper(structure_img_smooth, s2_param_bright)
bw_dark = dot_2d_slice_by_slice_wrapper(1-structure_img_smooth, s2_param_dark)

bw_merge = np.logical_or(bw, bw_extra)
bw_merge[bw_dark>0]=0


minArea = 5 * (1000 ** 2) / (pixel_size_in_nano ** 2) * (downsample_percentage ** 2) / (100 ** 2)

seg = remove_small_objects(bw_merge>0, min_size=minArea, connectivity=1, in_place=False)
seg = np.squeeze(seg, axis=0)
final_seg = label(seg * 1, background=0, connectivity=2)

imsave(join(file_dir, 'cell_mask_ACSS_classic.png'), final_seg)
pickle.dump(final_seg, bz2.BZ2File(join(file_dir, 'cell_mask_ACSS_classic.pkl'), 'w'))



