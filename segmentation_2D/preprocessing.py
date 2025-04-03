from pathlib import Path
from skimage.io import imread, imsave
import tifffile
import os
import numpy as np
import xml.etree.ElementTree as ET
#from skimage.measure import block_reduce
"""
FUNCTIONS TO READ IMAGE AND METADATA 
Author: Haoran Chen
Version: 1.3 February 14, 2025 R.F.Murphy
Version: 1.5.1 March 31m 2025 R.F.Murphy
         Replace -1 in crop limits with full dimension
"""

def extract_voxel_size_from_tiff(file_path):

    try:
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
    except:
        physical_size_x = 0.2
        physical_size_y = 0.2
        physical_size_z = 1.0

    return (physical_size_x, physical_size_y, physical_size_z)

def get_channel_names(img_dir, channel_names=None):

    if channel_names:
        channel_names = open(channel_names, 'r').read().splitlines()
        print(f"Channel names read from file: {channel_names}")
        return channel_names
    try:
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
        print(f"Channel names read from file: {channel_names}")
    except:
        channel_names = ["Ch1","Ch2"]
        print(f"Channel names assumed: {channel_names}")
    if not channel_names:
        channel_names = ["Ch1","Ch2"]

    return channel_names


def get_channel_intensity(marker_list, names, img):
    channel_intensity = np.zeros((img.shape[0], img.shape[2], img.shape[3]))
    for marker in marker_list:
        channel_idx = names.index(marker)
        channel_intensity = channel_intensity + img[:, channel_idx, :, :]
    return channel_intensity





def write_IMC_input_channels(img_file: Path, results_dir: Path, nucleus_channel_marker_list, cytoplasm_channel_marker_list,membrane_channel_marker_list,cl=None,channel_names=None):
    print(f"Crop limits: {cl}")
    image = imread(img_file)
    if cl != None:
        if(len(cl)!=6):
            print(f"Invalid crop limits: {cl}")
        else:
            #z, color/channel, y, x
            imageshap = image.shape
            #print(imageshap)
            if cl[1]<0:
                cl[1]=imageshap[0]
            if cl[3]<0:
                cl[3]=imageshap[2]
            if cl[5]<0:
                cl[5]=imageshap[3]
            image = image[cl[0]:cl[1],:,cl[2]:cl[3],cl[4]:cl[5]]
            print(f"Cropping image to shape {image.shape}")
    channel_names = get_channel_names(img_file, channel_names)
    
    nucleus_channel = get_channel_intensity(nucleus_channel_marker_list, channel_names, image)
    
    cytoplasm_channel = get_channel_intensity(cytoplasm_channel_marker_list, channel_names, image)
    
    membrane_channel = get_channel_intensity(membrane_channel_marker_list, channel_names, image)

    imsave(results_dir / 'nucleus.tif', nucleus_channel)
    imsave(results_dir / 'cytoplasm.tif', cytoplasm_channel)
    imsave(results_dir / 'membrane.tif', membrane_channel)

    #print(image.shape)
    
    return nucleus_channel, cytoplasm_channel, membrane_channel, image
