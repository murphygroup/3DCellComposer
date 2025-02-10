from pathlib import Path

from deepcell.applications import Mesmer
import numpy as np
from scipy.interpolate import interp1d

# Initialize TensorFlow
from tensorflow.compat.v1 import ConfigProto, InteractiveSession
from tensorflow.keras.models import load_model
import tensorflow as tf

"""
WRAPPER TO PERFORM 2D SEGMENTATIONS ALONG ALL AXES USING DEEPCELL
Author: Haoran Chen
Version: 1.1 December 14, 2023 R.F.Murphy
        Modify dimensions of the XZ and YZ to fully pad z with zeros  
"""

model_path = Path("/opt/.keras/models/0_12_9/MultiplexSegmentation")

def pad_to_multiple(arr, multiple=32):
	paddings = []
	for i in range(arr.ndim - 1):
		remainder = arr.shape[i] % multiple
		pad_amount = (multiple - remainder) % multiple
		paddings.append((0, pad_amount))
	paddings.append((0, 0))
	return np.pad(arr, paddings, mode='constant')

def interpolate_predictions(sampled_preds, total_slices, sampling_interval, interpolation_method='cubic', fill_value='extrapolate'):
    
    sampled_indices = np.arange(0, total_slices, sampling_interval)
    all_indices = np.arange(total_slices)
    
    # Options: 'cubic', 'previous', 'next', 'nearest'
    # Cubic often better for smooth transitions
    # interpolation_method = 'cubic'
    
    # Handle edge behavior
    bounds_error = False
    
	# options: 'extrapolate', 'constant'
    # constant often better for sharp transitions and 'extrapolate' often better for smooth transitions
    fill_value = 'extrapolate'
    
    pred_reshaped = sampled_preds.reshape(sampled_preds.shape[0], -1)
    interp_reshaped = np.zeros((total_slices, pred_reshaped.shape[1]))
    
    for i in range(pred_reshaped.shape[1]):
        f = interp1d(sampled_indices, 
                    pred_reshaped[:, i],
                    kind=interpolation_method,
                    bounds_error=bounds_error,
                    fill_value=fill_value)
        interp_reshaped[:, i] = f(all_indices)
    
    return interp_reshaped.reshape(total_slices, *sampled_preds.shape[1:])

def deepcell_segmentation_2D(im1, im2, axis, voxel_size, sampling_interval=3, chunk_size=100, interpolation_method='cubic', fill_value='extrapolate', dtype='float32'):
    """
    Perform 2D segmentations with slice sampling and interpolation
    """

    # z_slice_num = im1.shape[0]
    im1 = im1.astype(dtype)
    im2 = im2.astype(dtype)
    
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
    
    model = None
    if model_path.is_dir():
        model = load_model(model_path)
    
    config = ConfigProto()
    config.gpu_options.allow_growth = True
    session = InteractiveSession(config=config)
    config.gpu_options.per_process_gpu_memory_fraction = 0.9
    tf.compat.v1.keras.backend.set_session(tf.compat.v1.Session(config=config))
    app = Mesmer(model=model)
    
    print(f'Segmenting in {axis} direction with sampling interval {sampling_interval}...')
    
    # Process sampled slices
    sampled_slices = range(0, len(im), sampling_interval)
    predictions = []
    
    for chunk_start in range(0, len(sampled_slices), chunk_size):
        chunk_indices = list(sampled_slices)[chunk_start:chunk_start + chunk_size]
        chunk_predictions = []
        
        for i in chunk_indices:
            padded_slice = pad_to_multiple(im[i])
            pred = app.predict(np.expand_dims(padded_slice, 0), image_mpp=pixel_size, compartment='both')
            pred = pred[:, :im.shape[1], :im.shape[2], :]  # Crop back to original size
            chunk_predictions.append(pred[0])
        
        predictions.extend(chunk_predictions)
        print(f'Completed through sampled slice {chunk_indices[-1]} of {len(im)}')
    
    predictions = np.array(predictions)
    
    # Interpolate between sampled slices
    print('Interpolating between slices...')
    labeled_image = interpolate_predictions(predictions, len(im), sampling_interval, interpolation_method, fill_value)
    
    # Extract and rotate masks
    cell_mask = labeled_image[..., 0]
    nuc_mask = labeled_image[..., 1]
    
    if axis == 'XZ':
        cell_mask = np.rot90(cell_mask, k=-1, axes=(2, 0))
        nuc_mask = np.rot90(nuc_mask, k=-1, axes=(2, 0))
    elif axis == 'YZ':
        cell_mask = np.rot90(cell_mask, k=-1, axes=(1, 0))
        nuc_mask = np.rot90(nuc_mask, k=-1, axes=(1, 0))
    
    return cell_mask, nuc_mask

# def deepcell_segmentation_2D(im1, im2, axis, voxel_size):
# 	z_slice_num = im1.shape[0]
	
# 	if axis == 'XY':
# 		pixel_size = float(voxel_size[0])
# 	elif axis == 'XZ':
# 		im1 = np.rot90(im1, k=1, axes=(2, 0))
# 		im2 = np.rot90(im2, k=1, axes=(2, 0))
# 		pixel_size = float(voxel_size[2])
# 	elif axis == 'YZ':
# 		im1 = np.rot90(im1, k=1, axes=(1, 0))
# 		im2 = np.rot90(im2, k=1, axes=(1, 0))
# 		pixel_size = float(voxel_size[2])
	
# 	im = np.stack((im1, im2), axis=-1)
	
# 	if axis == 'XY':
# 		pass
# 	elif axis == 'XZ':
# 		# v1.0 300 was subtracted for XZ and YZ padding for faster running
# 		# v1.1 fully padded z for generalizability
# 		im_zeros = np.zeros((im.shape[0], im.shape[1], im.shape[0]-im.shape[2], im.shape[3])) # patch for DL network
# 		im = np.dstack((im, im_zeros))
# 	elif axis == 'YZ':
# 		im_zeros = np.zeros((im.shape[0], im.shape[2]-im.shape[1], im.shape[2], im.shape[3]))
# 		im = np.hstack((im, im_zeros))
# 	#print(im.shape)
	
	
# 	from tensorflow.compat.v1 import ConfigProto
# 	from tensorflow.compat.v1 import InteractiveSession
# 	from tensorflow.keras.models import load_model
# 	import tensorflow as tf

# 	model = None
# 	if model_path.is_dir():
# 		model = load_model(model_path)

# 	config = ConfigProto()
# 	config.gpu_options.allow_growth = True
# 	session = InteractiveSession(config=config)
# 	config.gpu_options.per_process_gpu_memory_fraction = 0.9
# 	tf.compat.v1.keras.backend.set_session(tf.compat.v1.Session(config=config))
# 	app = Mesmer(model=model)

# 	print('Segmenting in ',axis,' direction...')
# 	istep=max(1,int(len(im)/10))
# 	for i in range(len(im)):
# 		#print(tf.__version__)
# 		if i == 0:
# 			labeled_image = app.predict(np.expand_dims(im[i], 0), image_mpp=pixel_size, compartment='both')
# 		else:
# 			labeled_image = np.vstack((labeled_image, app.predict(np.expand_dims(im[i], 0), image_mpp=pixel_size, compartment='both')))
			
# 		if i>0 and (i%istep==0 or i==len(im)-1):
# 			print(f'Completed through slice {i}')
# 	if axis == 'XY':
# 		cell_mask = labeled_image[:, :, :, 0]
# 		nuc_mask = labeled_image[:, :, :, 1]
# 	elif axis == 'XZ':
# 		cell_mask = labeled_image[:, :, :z_slice_num, 0]
# 		nuc_mask = labeled_image[:, :, :z_slice_num, 1]
# 	elif axis == 'YZ':
# 		cell_mask = labeled_image[:, :z_slice_num, :, 0]
# 		nuc_mask = labeled_image[:, :z_slice_num, :, 1]
		
# 	return cell_mask, nuc_mask
