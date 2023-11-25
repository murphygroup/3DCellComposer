from deepcell.applications import Mesmer
import numpy as np


def deepcell_segmentation_2D(im1, im2, axis):
	z_slice_num = im1.shape[0]
	
	if axis == 'XY':
		pass
		pixel_size = 1
	elif axis == 'XZ':
		im1 = np.rot90(im1, k=1, axes=(2, 0))
		im2 = np.rot90(im2, k=1, axes=(2, 0))
		pixel_size = 2
	elif axis == 'YZ':
		im1 = np.rot90(im1, k=1, axes=(1, 0))
		im2 = np.rot90(im2, k=1, axes=(1, 0))
		pixel_size = 2
	
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
	import tensorflow as tf
	#
	config = ConfigProto()
	config.gpu_options.allow_growth = True
	session = InteractiveSession(config=config)
	config.gpu_options.per_process_gpu_memory_fraction = 0.9
	tf.compat.v1.keras.backend.set_session(tf.compat.v1.Session(config=config))
	app = Mesmer()
	
	for i in range(len(im)):
		print(f'Segmenting slice {i}')
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
		
	return cell_mask, nuc_mask

