# 3DCellComposer - A Versatile Pipeline Utilizing 2D Cell Segmentation Methods for 3D Cell Segmentation
Haoran Chen, Ted Chang, Matthew Ruffalo and Robert F. Murphy\
Carnegie Mellon University\
V1.5.5 July 8, 2025

3DCellComposer is a versatile, open-source software designed as a general solution for 3D cell segmentation. It allows users to choose any existing 2D segmentation model appropriate for their tissue or cell type(s) without requiring any additional training. Moreover, we have enhanced our CellSegmentationEvaluator quality evaluation tool to support 3D images. It allows users to compare and select the most suitable 2D segmentation models for 3D tasks, without the need for human annotations to assess performance.

It is available as a full-featured GitHub repository from <https://github.com/murphygroup/3DCellComposer>, and as this PyPI python package that provides a simplified implementation that uses just DeepCell as the 2D segmenter.  **Note that this package is significantly slower than (and does not support some of the options of) the version available from GitHub.**

Reference: Haoran Chen and Robert F. Murphy (2025) 3DCellComposer - A Versatile Pipeline Utilizing 2D Cell Segmentation Methods for 3D Cell Segmentation. Under review.

## Using the PyPI package

To use the package, 
```bash
pip install ThreeDCellComposer
```

To call the function in python,
```bash
from ThreeDCellComposer.ThreeDCellComposer import ThreeDCellComposer

ThreeDCellComposer(image_path,nucleus_channel_marker_list,cytoplasm_channel_marker_list,membrane_channel_marker_list,segmentation_method)
```
where the channel lists are strings consisting of the names (not numbers) of the channels to be used for segmentation. Only 'deepcell' is supported as the segmentation_method by the PyPI package at this time. Additional optional arguments are described below.

'run_3DCellComposerUsingPackage.py' is an example python script that calls ThreeDCellposer. 


## Contact

Robert F. Murphy - murphy@cmu.edu\
Haoran Chen - haoran.chen@stjude.org

