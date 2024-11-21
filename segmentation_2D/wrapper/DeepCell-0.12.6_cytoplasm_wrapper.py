from skimage.io import imread
from skimage.io import imsave
from skimage.color import rgb2gray
from deepcell.applications import Mesmer
# from deepcell.applications import NuclearSegmentation
# from deepcell.applications import CytoplasmSegmentation
import os
from pathlib import Path
from os.path import join
import numpy as np
import sys
import bz2
import pickle
import deepcell

model_path = Path("/opt/.keras/models/0_12_9/MultiplexSegmentation")

file_dir = sys.argv[1]
axis = sys.argv[2]
voxel_size = sys.argv[3]
im1 = imread(join(file_dir, 'nucleus.tif'))
im2 = imread(join(file_dir, 'cytoplasm.tif'))
z_slice_num = im1.shape[0]

if axis == 'XY':
	pixel_size = float(voxel_size[0])
elif axis == 'XZ':
	im1 = np.rot90(im1, k=1, axes=(2, 0))
	im2 = np.rot90(im2, k=1, axes=(2, 0))
	pixel_size = float(voxel_size[2])
elif axis == 'YZ':
	im1 = np.rot90(im1, k=1, axes=(1, 0))
	im2 = np.rot90(im2, k=1, axes=(1, 0))
	pixel_size = float(voxel_size[2])


im = np.stack((im1, im2), axis=-1)


if axis == 'XY':
	pass
elif axis == 'XZ':
	im_zeros = np.zeros((im.shape[0], im.shape[1], im.shape[0]-im.shape[2]-300, im.shape[3]))
	im = np.dstack((im, im_zeros))
elif axis == 'YZ':
	im_zeros = np.zeros((im.shape[0], im.shape[2]-im.shape[1]-300, im.shape[2], im.shape[3]))
	im = np.hstack((im, im_zeros))


from tensorflow.compat.v1 import ConfigProto
from tensorflow.compat.v1 import InteractiveSession
from tensorflow.keras.models import load_model
import tensorflow as tf

model = None
if model_path.is_dir():
	model = load_model(model_path)

config = ConfigProto()
config.gpu_options.allow_growth = True
session = InteractiveSession(config=config)
config.gpu_options.per_process_gpu_memory_fraction = 0.9
tf.compat.v1.keras.backend.set_session(tf.compat.v1.Session(config=config))
app = Mesmer(model=model)

for i in range(len(im)):
	if i == 0:
		labeled_image = app.predict(np.expand_dims(im[i], 0), image_mpp=pixel_size, compartment='both')
	else:
		labeled_image = np.vstack((labeled_image, app.predict(np.expand_dims(im[i], 0), image_mpp=pixel_size, compartment='both')))
		
if axis == 'XY':
	cell_mask = labeled_image[:, :, :, 0]
	nuc_mask = labeled_image[:, :, :, 1]
elif axis == 'XZ':
	cell_mask = labeled_image[:, :, :z_slice_num, 0]
	nuc_mask = labeled_image[:, :, :z_slice_num, 1]
elif axis == 'YZ':
	cell_mask = labeled_image[:, :z_slice_num, :, 0]
	nuc_mask = labeled_image[:, :z_slice_num, :, 1]

mask_dir = bz2.BZ2File(join(file_dir, 'cell_mask_DeepCell-0.12.6_cytoplasm_' + axis + '.pkl'), 'wb')
pickle.dump(cell_mask, mask_dir)

nuc_mask_dir = bz2.BZ2File(join(file_dir, 'nuclear_mask_DeepCell-0.12.6_cytoplasm_' + axis + '.pkl'), 'wb')
pickle.dump(nuc_mask, nuc_mask_dir)
