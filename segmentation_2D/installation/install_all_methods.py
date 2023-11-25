import os
def install_segmentation_methods():
	dir_path = os.path.dirname(os.path.realpath(__file__))
	
	# print('installing DeepCell-0.12.6...')
	# os.system(f'bash {dir_path}/install_DeepCell-0.12.6.sh {dir_path}')
	#
	# print('installing Cellpose-2.2.2...')
	# os.system(f'bash {dir_path}/install_Cellpose-2.2.2.sh {dir_path}')
	
	# print('installing ACSS classic...')
	# os.system(f'bash {dir_path}/install_ACSS_classic.sh {dir_path}')
	
	print('installing CellProfiler...')
	os.system(f'bash {dir_path}/install_CellProfiler.sh {dir_path}')
	
	# print('installing CellSegm...')
	# os.system(f'bash {dir_path}/install_CellSegm.sh {os.path.dirname(dir_path)}/cellsegm')
	#
	# print('installing CellX...')
	# os.system(f'bash {dir_path}/install_CellX.sh {os.path.dirname(dir_path)}/CellX')
