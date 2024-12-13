from pathlib import Path

from deepcell.applications import Mesmer
import numpy as np

"""
WRAPPER TO PERFORM 2D SEGMENTATIONS ALONG ALL AXES USING DEEPCELL
Author: Haoran Chen
Version: 1.1 December 14, 2023 R.F.Murphy
        Modify dimensions of the XZ and YZ to fully pad z with zeros  
"""

model_path = Path("/opt/.keras/models/0_12_9/MultiplexSegmentation")

def deepcell_segmentation_2D(im1, im2, axis, voxel_size, subsample_rate=2):
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
	# 	# v1.0 300 was subtracted for XZ and YZ padding for faster running
	# 	# v1.1 fully padded z for generalizability
		im_zeros = np.zeros((im.shape[0], im.shape[1], im.shape[0]-im.shape[2], im.shape[3])) # patch for DL network
		im = np.dstack((im, im_zeros))
		# im = np.concatenate([im] * int(abs(im.shape[0]-im.shape[2])/im.shape[2]), axis=2)  # This duplicates the Z-axis
	elif axis == 'YZ':
		im_zeros = np.zeros((im.shape[0], im.shape[2]-im.shape[1], im.shape[2], im.shape[3]))
		im = np.hstack((im, im_zeros))
		# im = np.concatenate([im] * int(abs(im.shape[2]-im.shape[1])/im.shape[2]), axis=1)  # This duplicates the Z-axis
	# print(im.shape)
	
	
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

	print('Segmenting in ',axis,' direction...')
	istep=max(1,int(len(im)/10))
	print(im.shape)

	batch_size = 200
	predictions = []

	#subsample 
	# subsample_im = im[::subsample_rate]

	#if len(im) < batch size
	if len(im) < batch_size:
		batch_size = len(im)
 
	for i in range(0, len(im), batch_size):
		batch = im[i:i + batch_size]
		batch_preds = app.predict(batch, batch_size=batch_size, image_mpp=pixel_size, compartment='both')
		predictions.append(batch_preds)

		if i > 0 and (i%istep==0 or i==len(im) - 1):
			print(f'Completed through slice {i}')

	# Concatenate all predictions along the first dimension
	labeled_image = np.concatenate(predictions, axis=0)

	# for i in range(0, len(im), batch_size):
	# 	#print(tf.__version__)
	# 	# print(np.amax(im[i]))
	# 	if i == 0:
	# 		labeled_image = app.predict(np.expand_dims(im[i], 0), image_mpp=pixel_size, compartment='both')
	# 	else:
	# 		labeled_image = np.vstack((labeled_image, app.predict(np.expand_dims(im[i], 0), image_mpp=pixel_size, compartment='both')))
	# 	print(np.amax(labeled_image))
			
	# 	if i>0 and (i%istep==0 or i==len(im)-1):
	# 		print(f'Completed through slice {i}')
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
