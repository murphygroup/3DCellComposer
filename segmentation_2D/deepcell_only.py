from pathlib import Path

from deepcell.applications import Mesmer
import numpy as np
import math
import pickle
import os
from datetime import datetime
import time

# Initialize TensorFlow
from tensorflow.compat.v1 import ConfigProto, InteractiveSession
from tensorflow.keras.models import load_model
import tensorflow as tf

"""
WRAPPER TO PERFORM 2D SEGMENTATIONS ALONG ALL AXES USING DEEPCELL
Author: Haoran Chen
Version: 1.1 December 14, 2023 R.F.Murphy
        Modify dimensions of the XZ and YZ to fully pad z with zeros
Version: 1.3 February 14, 2025 R.F.Murphy
        copy segmented slices to fill between them
        save prevously segmented slices for efficiency
Version: 1.5 March 7, 2025 R.F.Murphy
        allow adjustable minimum slice padding
"""

model_path = Path("/opt/.keras/models/0_12_9/MultiplexSegmentation")

def pad_to_multiple(arr, multiple=512):
	paddings = []
	for i in range(arr.ndim - 1):
		pad_amount = 0
		if arr.shape[i] < multiple:
			remainder = arr.shape[i] % multiple
			pad_amount = (multiple - remainder) % multiple
		paddings.append((0, pad_amount))
	paddings.append((0, 0))
	return np.pad(arr, paddings, mode='constant')

def fill_in_slices(sampled_preds, total_slices):
    siz=sampled_preds.shape
    #print(siz)
    nrepeat = math.ceil(total_slices / siz[0])
    #print(nrepeat)
    full = np.repeat(sampled_preds,nrepeat,axis=0)
    #print(full.shape)
    return full[0:total_slices,:,:]

def fill_in_slices3D(sampled_preds, desired_shape):
    siz=sampled_preds.shape
    print(siz)
    print(desired_shape)
    nrepeat0 = math.ceil(desired_shape[0] / siz[0])
    print(nrepeat0)
    next0 = np.repeat(sampled_preds,nrepeat0,axis=0)
    print(next0.shape)
    nrepeat1 = math.ceil(desired_shape[1] / siz[1])
    print(nrepeat1)
    next1 = np.repeat(next0,nrepeat1,axis=1)
    print(next1.shape)
    nrepeat2 = math.ceil(desired_shape[2] / siz[2])
    print(nrepeat2)
    full = np.repeat(next1,nrepeat2,axis=2)
    print(full.shape)
    return full[0:desired_shape[0],0:desired_shape[1],0:desired_shape[2]]

def deepcell_segmentation_2D(im1, im2, axis, voxel_size, sampling_interval=3, chunk_size=100, results_path=Path('results'), maxima_threshold=0.075, interior_threshold=0.2, compartment='both', min_slice_padding='512', dtype='float32'):
    """
    Perform 2D segmentations with slice sampling
    """

    #print(type(im1),type(im2))
    #print(np.min(im1),np.max(im1))
    #print(np.min(im2),np.max(im2))
    #dtype = 'int16'
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

    ssfilename = results_path / ('saved_segmentations' + axis + '.pkl')
    if os.path.exists(ssfilename):
        with open(ssfilename, 'rb') as ss_file:
            saved_segmentations = pickle.load(ss_file)
    else:
        saved_segmentations = [None] * len(im)

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
    desired_slices = range(0, len(im), sampling_interval)
    #print(desired_slices)
    predictions = []

    for i in desired_slices:
        if not isinstance(saved_segmentations[i],np.ndarray):
            padded_slice = pad_to_multiple(im[i],min_slice_padding)
            if i==desired_slices[0]:
                print(f"Padded shape: {padded_slice.shape}")
            tstart = time.time()
            #pred = app.predict(np.expand_dims(padded_slice, 0), image_mpp=pixel_size, compartment='both')
            #pred = app.predict(np.expand_dims(padded_slice, 0), image_mpp=pixel_size, compartment='both', postprocess_kwargs_whole_cell={'maxima_threshold': maxima_threshold})
            #expslice = np.expand_dims(padded_slice, 0)
            #print(f"Expanded slice shape: {expslice.shape}")
            #pred = app.predict(np.expand_dims(padded_slice, 0), image_mpp=pixel_size, compartment='whole-cell', postprocess_kwargs_whole_cell={'maxima_threshold': 0.2})
            pred = app.predict(np.expand_dims(padded_slice, 0), image_mpp=pixel_size, compartment=compartment, postprocess_kwargs_whole_cell={'interior_threshold': interior_threshold, 'maxima_threshold': maxima_threshold})
            #setting to whole-cell only returns x,y,1 instead of x,y,2 so dup
            #print(pred.shape)
            pred = np.repeat(pred,2,axis=3)
            #print(pred.shape)
            #pred = app.predict(np.expand_dims(padded_slice, 0), image_mpp=pixel_size, compartment='both', postprocess_kwargs_whole_cell={'interior_threshold': 0.5})
            #pred = app.predict(np.expand_dims(padded_slice, 0), image_mpp=pixel_size, compartment='cytoplasm', postprocess_kwargs_whole_cell={'interior_threshold': 0.5})
            pred = pred[:, :im.shape[1], :im.shape[2], :]  # Crop back to orig
            if not i%chunk_size:
                print(f'Segmented slice {i} of {len(im)}')
                #print(f"max in slice {i} is {np.max(pred)}")
                #print(pred.shape)
            if i == desired_slices[0]:
                tend = time.time()
                etime = tend-tstart
                esttime = etime*(len(desired_slices)-1)/60.
                print(f"Elapsed time: {round(etime/60.,2)} min; estimated remaining time for this axis: {round(esttime,2)} min")
            saved_segmentations[i]=pred[0]
        predictions.append(saved_segmentations[i])
    #print(len(predictions),len(saved_segmentations))

    with open(ssfilename, 'wb') as ss_file:
        pickle.dump(saved_segmentations, ss_file)

    predictions = np.array(predictions)
    #print(predictions.shape)
    
    if sampling_interval > 1:
        print('Filling between sampled slices...')
    #print(predictions.shape)
    labeled_image = fill_in_slices(predictions, len(im))
    #print(labeled_image.shape)
    
    # Extract and rotate masks
    cell_mask = labeled_image[..., 0]
    nuc_mask = labeled_image[..., 1]
    #cell_mask = predictions[..., 0]
    #nuc_mask = predictions[..., 1]
    
    if axis == 'XZ':
        cell_mask = np.rot90(cell_mask, k=-1, axes=(2, 0))
        nuc_mask = np.rot90(nuc_mask, k=-1, axes=(2, 0))
    elif axis == 'YZ':
        cell_mask = np.rot90(cell_mask, k=-1, axes=(1, 0))
        nuc_mask = np.rot90(nuc_mask, k=-1, axes=(1, 0))
    
    return cell_mask, nuc_mask

