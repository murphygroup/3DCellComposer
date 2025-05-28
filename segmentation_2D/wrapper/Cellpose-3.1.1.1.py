import numpy as np
from os.path import join
from cellpose.io import imread, imsave
#import bz2
#import pickle
import sys
from cellpose import models

file_dir = sys.argv[1]

im1 = imread(join(file_dir, 'cytoplasm.tif'))
im2 = imread(join(file_dir, 'nucleus.tif'))
#im = np.stack((im1, im2))
im = np.stack((im1, im2),axis=2)

###model = models.Cellpose(gpu=True)

#call sequence for default model
###masks = model.eval(im, flow_threshold=0.4, cellprob_threshold=0.0)
###masks = masks[0]
#call sequence for 3.1.1.1
masks, flows, styles = models.CellposeModel(model_type='tissuenet_cp3',gpu=True).eval(im, diameter=0, channels=[0,0])
#call sequence for 4.0.4
#masks, flows, styles = models.CellposeModel(gpu=True).eval(im)
imsave(join(file_dir, 'mask_Cellpose-3.1.1.1.tif'), masks)
#mask_dir = bz2.BZ2File(join(file_dir, 'mask_Cellpose-3.1.1.1.pickle'), 'wb')
#pickle.dump(masks, mask_dir)

###masksnuc = model.eval(im2, flow_threshold=0.4, cellprob_threshold=0.0)
###masksnuc = masksnuc[0]
masksnuc, flows, styles = models.CellposeModel(model_type='tissuenet_cp3').eval(im2, diameter=0)
#masksnuc, flows, styles = models.CellposeModel(gpu=True).eval(im2)
imsave(join(file_dir, 'nuclear_mask_Cellpose-3.1.1.1.tif'), masksnuc)
#nuc_mask_dir = bz2.BZ2File(join(file_dir, 'nuclear_mask_Cellpose-3.1.1.1.pickle'), 'wb')
#pickle.dump(masksnuc, nuc_mask_dir)
