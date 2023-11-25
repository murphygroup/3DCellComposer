from skimage.io import imread
import sys
import numpy as np
from os.path import join
import bz2
import pickle

from scipy.ndimage import label

def convert_to_indexed(img):
	unique_cell = sorted(np.unique(img))
	print(len(unique_cell))
	n_cell = len(np.unique(img))
	for i in range(1, n_cell):
		img[np.where(img == unique_cell[i])] = i
	return img

if __name__ == '__main__':
	file_dir = sys.argv[1]
	try:
		img = np.load(join(file_dir, 'cell_mask_' + sys.argv[2] + '.npy'))
	except:
		img = imread(join(file_dir, 'cell_mask_' + sys.argv[2] + '.png'), as_gray=True)


	if sys.argv[2] == 'CellX' or sys.argv[2] == 'CellSegm':
		img = label(img)[0]

	else:
		img = convert_to_indexed(img)

	save_dir = bz2.BZ2File(join(file_dir, 'cell_mask_' + sys.argv[2] + '.pkl'), 'w')
	pickle.dump(img.astype(int), save_dir)
	
