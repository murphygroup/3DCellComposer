from skimage.io import imread, imsave
import tifffile
import os
import numpy as np
import xml.etree.ElementTree as ET

'''
version 1.3: February 19, 2025 R.F.Murphy
    add downsample support

def extract_voxel_size_from_tiff(file_path):
    # Read OME-TIFF metadata
    with tifffile.TiffFile(file_path) as tif:
        metadata = tif.ome_metadata

    # Parse the XML metadata
    root = ET.fromstring(metadata)

    # Initialize sizes
    physical_size_x = physical_size_y = physical_size_z = None

    # Iterate over all elements to find the first instance with physical sizes
    for elem in root.iter():
        if 'PhysicalSizeX' in elem.attrib:
            physical_size_x = elem.get('PhysicalSizeX')
        if 'PhysicalSizeY' in elem.attrib:
            physical_size_y = elem.get('PhysicalSizeY')
        if 'PhysicalSizeZ' in elem.attrib:
            physical_size_z = elem.get('PhysicalSizeZ')
        if physical_size_x and physical_size_y and physical_size_z:
            break

    return (physical_size_x, physical_size_y, physical_size_z)

def get_channel_names(img_dir):
    with tifffile.TiffFile(img_dir) as tif:
        # Get the ImageDescription which contains the XML metadata
        description = tif.pages[0].tags['ImageDescription'].value

    # Parse the XML metadata
    root = ET.fromstring(description)

    # Find all Channel elements and extract their names
    channel_names = []
    for elem in root.iter():
        if 'Channel' in elem.tag:
            name = elem.get('Name')
            if name:
                channel_names.append(name)

    return channel_names


def get_channel_intensity(marker_list, names, img):
    channel_intensity = np.zeros((img.shape[0], img.shape[2], img.shape[3]))
    for marker in marker_list:
        channel_idx = names.index(marker)
        channel_intensity = channel_intensity + img[:, channel_idx, :, :]
    return channel_intensity





def write_IMC_input_channels(img_dir, nucleus_channel_marker_list, cytoplasm_channel_marker_list,
                             membrane_channel_marker_list,downsample_vector):
    image = imread(img_dir)
    channel_names = get_channel_names(img_dir)
    
    # nucleus_channel_marker_list = ['Ir191']
    nucleus_channel = get_channel_intensity(nucleus_channel_marker_list, channel_names, image)
    
    # cytoplasm_channel_marker_list = ['In115', 'Y89', 'Tb159']
    cytoplasm_channel = get_channel_intensity(cytoplasm_channel_marker_list, channel_names, image)
    
    # membrane_channel_marker_list = ['La139', 'Pr141', 'Eu151', 'Gd160', 'Dy162']
    membrane_channel = get_channel_intensity(membrane_channel_marker_list, channel_names, image)

    nucleus_channel = skimage.measure.block_reduce(nucleus_channel, block_size=downsample_vector,func=np.max)
    cytoplasm_channel = skimage.measure.block_reduce(cytoplasm_channel, block_size=downsample_vector,func=np.max)
    membrane_channel = skimage.measure.block_reduce(membrane_channel, block_size=downsample_vector,func=np.max)
    
    imsave(f'{os.path.dirname(img_dir)}/nucleus.tif', nucleus_channel)
    imsave(f'{os.path.dirname(img_dir)}/cytoplasm.tif', cytoplasm_channel)
    imsave(f'{os.path.dirname(img_dir)}/membrane.tif', membrane_channel)
    
    return nucleus_channel, cytoplasm_channel, membrane_channel, image
