import numpy as np
from CellSegmentationEvaluator.CSE3D import CSE3D
from CellSegmentationEvaluator.functions import thresholding

"""
FUNCTION USED BY 3DCELLCOMPOSER TO CALCULATE SEGMENTATION EVALUATION 
STATISTICS FOR A SINGLE 3D IMAGE AND CELL AND NUCLEAR MASKS
Author: Haoran Chen
Version: 1.1 December 14, 2023 R.F.Murphy, Haoran Chen
         Modify KMeans section to avoid doing KMeans with n_clusters=1
         Remove saving intermediate binary image file to avoid error in the next run
         1.2 December 24, 2023 R.F.Murphy
         Use CSE3D from CellSegmentationEvaluator
         1.3 March 15, 2025 R.F.Murphy
         remove unnecessary print statements
"""

def seg_evaluation_3D(cell_matched_mask,
                      nuclear_matched_mask,
                      nucleus,
                      cytoplasm,
                      membrane,
                      img_channels,
                      voxel_size,
                      pca_dir):

	
	cell_outside_nucleus_mask = cell_matched_mask - nuclear_matched_mask
	metric_mask = np.expand_dims(cell_matched_mask, 0)
	metric_mask = np.vstack((metric_mask, np.expand_dims(nuclear_matched_mask, 0)))
	metric_mask = np.vstack((metric_mask, np.expand_dims(cell_outside_nucleus_mask, 0)))

	#print('in single_method_eval')
	#print(np.amax(nucleus),np.amax(cytoplasm),np.amax(membrane))
	img_input_channel = nucleus + cytoplasm + membrane
	#print(img_input_channel.shape,np.amax(img_input_channel))
	#for slice in range(img_input_channel.shape[0]):
		#print(np.amin(img_input_channel[slice,:,:]),np.amax(img_input_channel[slice,:,:]),np.sum(img_input_channel[slice,:,:]))
	img4thresh = np.stack([thresholding(img_input_channel[slice,:,:]) for slice in range(img_input_channel.shape[0])])
	#img4thresh = np.sign(img_img4thresh)
	#print('min,max,sum for img4thresh')
	#print(type(img4thresh))
	#print(np.amin(img4thresh),np.amax(img4thresh),np.sum(img4thresh))

	# 3D IMC has dim order ZCYX
	img_channels = np.transpose(img_channels, (1, 0, 2, 3)) #TODO: automatically detect dim order
	#print('Channel image dimensions=',img_channels.shape)

	PCA_model = "3Dv1.6"

	vox_size = float(voxel_size[0]) * float(voxel_size[1]) * float(voxel_size[2])

	#metrics = CSE3D(img_channels, metric_mask, PCA_model, img4thresh, vox_size, (1,2,5), (20000,1000))
	metrics = CSE3D(img_channels, metric_mask, PCA_model, img4thresh, vox_size, (1,2,5), (20000,0))
	weighted_score = metrics["QualityScore"]
        
	return weighted_score, metrics
