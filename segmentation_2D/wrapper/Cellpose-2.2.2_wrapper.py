import numpy as np
import time, os, sys
from os.path import join
from skimage.io import imread
from skimage.io import imsave
import bz2
import pickle
import sys

use_GPU = 1
print('GPU activated? %d'%use_GPU)

file_dir = sys.argv[1]
axis = sys.argv[2]


im1 = imread(join(file_dir, 'cytoplasm.tif'))
im2 = imread(join(file_dir, 'nucleus.tif'))
if axis == 'XY':
	pass
	pixel_size = 0.1083
elif axis == 'XZ':
	im1 = np.rot90(im1, k=1, axes=(2, 0))
	im2 = np.rot90(im2, k=1, axes=(2, 0))
	pixel_size = 0.29
elif axis == 'YZ':
	im1 = np.rot90(im1, k=1, axes=(1, 0))
	im2 = np.rot90(im2, k=1, axes=(1, 0))
	pixel_size = 0.29
	
pixel_size_ratio = pixel_size / 1
im = np.stack((im1, im2))
from cellpose import models

# DEFINE CELLPOSE MODEL
model = models.CellposeModel(gpu=use_GPU, model_type='cyto')
channels = [0, 1]


import torch
torch.cuda.empty_cache()
masks_3D_pieces = list()
for z in range(im.shape[1]):
	masks, flows, styles = model.eval(im[:,z,:,:], diameter=10/pixel_size_ratio, channels=channels)
	masks_3D_pieces.append(masks)
masks_3D = np.stack(masks_3D_pieces, axis=0)
pickle.dump(masks_3D, bz2.BZ2File(f'{file_dir}/cell_mask_Cellpose-2.2.2_{axis}.pkl','w'))

model_nuc = models.CellposeModel(gpu=use_GPU, model_type='cyto')
channels = [0, 0]
im2 = np.expand_dims(im2, axis=0)
nuclear_masks_3D_pieces = list()
for z in range(im.shape[1]):
	masks, flows, styles = model.eval(im2[:,z,:,:], diameter=8/pixel_size_ratio, channels=channels)
	nuclear_masks_3D_pieces.append(masks)

nuclear_masks_3D = np.stack(nuclear_masks_3D_pieces, axis=0)

pickle.dump(nuclear_masks_3D, bz2.BZ2File(f'{file_dir}/nuclear_mask_Cellpose-2.2.2_{axis}.pkl','w'))

print('segmentation done')
