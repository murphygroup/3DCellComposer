import numpy as np
import glob
from os.path import join
import os
from skimage.io import imread
import bz2
import pickle

def extract_number(s):
	return int(''.join(filter(str.isdigit, s)))

def combine_slices():
	data_dir = '/data/3D/IMC_3D/florida-3d-imc/'
	img_3D_dir_list = sorted(glob.glob(join(data_dir, '**', 'random_gaussian_*'), recursive=True)) + sorted(glob.glob(join(data_dir, '**', 'original'), recursive=True))
	axes = ['XY', 'XZ', 'YZ']
	# methods = ['aics_classic', 'cellpose', 'CellProfiler', 'CellX', 'cellsegm', 'artificial']
	# methods = ['CellProfiler', 'aics_classic', 'CellX', 'cellsegm', 'artificial']
	# methods = ['CellProfiler', 'aics_classic', 'CellX', 'cellsegm']
	methods = ['artificial']
	for img_dir in img_3D_dir_list:
		print(img_dir)
		img_shape = imread(join(img_dir, 'nucleus.tif')).shape

		for method in methods:
			print(method)
			for axis in axes:
				if axis == 'XY':
					slice_shape = (img_shape[1], img_shape[2])
				elif axis == 'YZ':
					slice_shape = (img_shape[0], img_shape[2])
				elif axis == 'XZ':
					slice_shape = (img_shape[1], img_shape[0])
				# if not os.path.exists(f'{img_dir}/mask_{method}_{axis}.pkl'):
				if True:
					img_3D_pieces = list()
					slice_dir_list = sorted(glob.glob(f'{img_dir}/slices/slice_{axis}*', recursive=True), key=extract_number)
					for slice_dir in slice_dir_list:
						try:
							img_slice = pickle.load(bz2.BZ2File(f'{slice_dir}/mask_{method}.pkl', 'r'))
						except:
							print(slice_dir)
							img_slice = np.zeros(slice_shape)
						img_3D_pieces.append(img_slice)
					img_3D = np.stack(img_3D_pieces, axis=0)
					pickle.dump(img_3D, bz2.BZ2File(f'{img_dir}/mask_{method}_{axis}.pkl', 'w'))
				else:
					print(f'{img_dir}/mask_{method}_{axis}.pkl')
					
def split_slices():
	data_dir = '/data/3D/IMC_3D/florida-3d-imc/'
	img_3D_dir_list = sorted(glob.glob(join(data_dir, '**', 'original'), recursive=True)) + sorted(glob.glob(join(data_dir, '**', 'random_gaussian_*'), recursive=True))
	img_names = ['nucleus', 'cytoplasm', 'membrane']
	axes = ['XY', 'XZ', 'YZ']
	for img_dir in img_3D_dir_list:
		for img_name in img_names:
			print(join(img_dir, img_name + '.tif'))
			img = imread(join(img_dir, img_name + '.tif'))
			for axis in axes:
				if axis == 'XY':
					img_rotated = img.copy()
				elif axis == 'XZ':
					img_rotated = np.rot90(img, k=1, axes=(2, 0)).copy()
				elif axis == 'YZ':
					img_rotated = np.rot90(img, k=1, axes=(1, 0)).copy()
				print(axis)
				print(img_rotated.shape)
				for slice_index in range(img_rotated.shape[0]):
					img_slice = img_rotated[slice_index]
					os.system('rm -rf ' + join(img_dir, str(axis) + '_' + str(slice_index)))
					slice_dir = join(img_dir, 'slices', 'slice_' + str(axis) + '_' + str(slice_index))
					if not os.path.exists(slice_dir):
						os.makedirs(slice_dir)
					save_dir = join(slice_dir, img_name + '.tif')
					
					
def segmentation_single_method(method_name):
