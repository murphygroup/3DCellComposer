import numpy as np
import glob
from os.path import join
import os
from skimage.io import imread, imsave
import bz2
import pickle
import subprocess


def extract_number(s):
	return int(''.join(filter(str.isdigit, s)))

def combine_slices(img_dir, method):
	axes = ['XY', 'XZ', 'YZ']
	img_shape = imread(join(img_dir, 'nucleus.tif')).shape
	for axis in axes:
		if axis == 'XY':
			slice_shape = (img_shape[1], img_shape[2])
		elif axis == 'YZ':
			slice_shape = (img_shape[0], img_shape[2])
		elif axis == 'XZ':
			slice_shape = (img_shape[1], img_shape[0])
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
		

					
def split_slices(img_dir):
	img_names = ['nucleus', 'cytoplasm', 'membrane']
	axes = ['XY', 'XZ', 'YZ']
	for img_name in img_names:
		img = imread(join(img_dir, img_name + '.tif'))
		for axis in axes:
			if axis == 'XY':
				img_rotated = img.copy()
			elif axis == 'XZ':
				img_rotated = np.rot90(img, k=1, axes=(2, 0)).copy()
			elif axis == 'YZ':
				img_rotated = np.rot90(img, k=1, axes=(1, 0)).copy()
			for slice_index in range(img_rotated.shape[0]):
				img_slice = img_rotated[slice_index]
				os.system('rm -rf ' + join(img_dir, str(axis) + '_' + str(slice_index)))
				slice_dir = join(img_dir, 'slices', 'slice_' + str(axis) + '_' + str(slice_index))
				if not os.path.exists(slice_dir):
					os.makedirs(slice_dir)
				save_dir = join(slice_dir, img_name + '.tif')
				imsave(save_dir, img_slice.astype(np.uint16))
					
def  segmentation_single_method(method, img_path):
	current_dir = os.getcwd()
	run_file = f'run_{method}.sh'
	if method in ['DeepCell-0.12.6_membrane', 'DeepCell-0.12.6_cytoplasm', 'Cellpose-2.2.2']:
		subprocess.run([f'{current_dir}/segmentation_2D/{run_file}', img_path])
	elif method in ['CellProfiler', 'ACSS_classic', 'CellX', 'CellSegm']:
		axes = ['XY', 'XZ', 'YZ']
		for axis in axes[:1]:
			slice_paths = sorted(glob.glob(f'{img_path}/slices/slice_{axis}*', recursive=True))
			if axis == 'XY':
				pixel_size = '1000'
			else:
				pixel_size = '2000'
			downsample_percentage = '100' # no downsample perturbation introduced
			for slice_path in slice_paths[:1]:
				print(slice_path, pixel_size)
				subprocess.run([f'{current_dir}/segmentation_2D/{run_file}', slice_path, downsample_percentage, pixel_size])

			

	