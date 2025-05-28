# 3DCellComposer - A Versatile Pipeline Utilizing 2D Cell Segmentation Methods for 3D Cell Segmentation
Haoran Chen, Ted Chang, Matthew Ruffalo and Robert F. Murphy\
Carnegie Mellon University\
V1.5.2 May 27, 2025

3DCellComposer is a versatile, open-source software designed as a general solution for 3D cell segmentation. It allows users to choose an existing 2D segmentation model appropriate for their tissue or cell type(s) without requiring any additional training. Moreover, we have enhanced our CellSegmentationEvaluator quality evaluation tool to support 3D images. It allows users to compare and select the most suitable 2D segmentation models for 3D tasks, without the need for human annotations to assess performance.

It is available as a full-featured GitHub repository, and as a python package on PyPi that provides an older, simplified implementation that uses just DeepCell as the 2D segmenter.

Reference: Haoran Chen and Robert F. Murphy (2025) 3DCellComposer - A Versatile Pipeline Utilizing 2D Cell Segmentation Methods for 3D Cell Segmentation. Under review.

Changes in V1.3: Options added to support subsampling of the slices along each axis.  Caching of slice segmentation results to avoid segmenting again.  Fixes to image padding when the number of z slices is too small.

Changes in V1.5: Added options for cropping initial image, block-wise voxel downsampling prior to segmentation,- setting specific DeepCell parameters (including whether to segment using cell channel, nuclear channel or both), and controlling the minimum padding

Changes in V1.5.2: Improved installation instructions and update environment.yml and requirements.txt.  Updated Cellpose support. Removed use of pint package.

## Using the repository version

### Overview
The `3DCellComposer` script processes multiplexed imaging data for cell segmentation. It requires specific inputs including the path to the image file and lists of markers for different cell components (nucleus, cytoplasm, and cell membrane). It currently can use DeepCell,  Cellpose or a user-provided method for 2D segmentation

### Installation

- Clone or download this GitHub repository.
- Install conda and pip, if not already installed.
- If planning to use DeepCell with 3DCellComposer (see note about required access token below), note that DeepCell requires 'tensorflow' which currently is not supported for python versions greater than 3.12.  If you are running a later version of python, you can create a conda environment that uses an earlier python version by doing
``bash
conda create --name 3DCellComposer python=3.9.13
conda activate 3DCellComposer
```
- 3DCellComposer can use either DeepCell or Cellpose to perform 2D segmentations.  Install the required packages one of two ways:
  (Recommended) Install required packages (note this adds them to your current environment):
``bash
pip install -r requirements.txt
```
      or
  (If your python version is 3.12 or below) Create a new conda environment from a specification file:
``bash
conda env create -f environment.yml
```
- The software was tested on Python 3.8 and 3.9 on Ubuntu 18.04.5 LTS and MacOS 14.7.4 and 14.7.6

- Note: Segmentation with DeepCell requires an access token that can be obtained at https://users.deepcell.org/.  The token needs to be added to your environment using ``bash
export DEEPCELL_ACCESS_TOKEN=<token_value>
```

### Command Structure
3DCellComposer is executed with a command like the following:

```bash
python run_3DCellComposer.py [image_path] [nucleus_markers] [cytoplasm_markers] [membrane_markers]
```
### Description of Required Arguments

1. **Image Path**: 
   - **Description**: Path to your multiplexed image file.
   - **Format**: String
   - **Example**: `/path/to/your/image.ome.tiff`

2. **Nucleus Channel Marker List**:
   - **Description**: A list of nuclear marker(s) used in the multiplexed image for segmentation.
   - **Format**: String of markers separated by commas, no spaces.
   - **Example**: For a single marker: `"Ir191"`, for multiple markers: `"In115,Y89,Tb159"`

3. **Cytoplasm Channel Marker List**:
   - **Description**: A list of cytoplasmic marker(s) used in the multiplexed image for segmentation.
   - **Format**: Similar to the nucleus channel marker list.
   - **Example**: `"La139,Pr141,Eu151"`

4. **Membrane Channel Marker List**:
   - **Description**: A list of cell membrane marker(s) used in the multiplexed image for segmentation.
   - **Format**: Similar to the nucleus and cytoplasm channel marker lists.
   - **Example**: `"Gd160,Dy162"`

### Most frequently used optional argument descriptions

**--skipYZ**
**Description**: Only do 2D segmentation on XY and XZ slices
**Format*: True or False
default=False

**--segmentation_method**
**Description**: Choose the 2D segmentation method.
      - "deepcell" - Employs DeepCell 2D segmentation, which achieved the highest performance in our evaluation.
      - "cellpose" - Employs Cellpose 2D segmentation.
      - "custom" - Utilizes a user-provided segmentation method, for which an empty wrapper with instructions is supplied in `3DCellComposer/segmentation_2D/wrapper/custom_segmentation_wrapper.py`.
**Format**: String (one of "deepcell", cellpose", "custom")
**Default**: "deepcell"

**--results_path**
**Description**: Specify path where results files will be written
**Format*: Path
default=Path('results')

**--crop_limits**
**Description**: Sets a region to crop the input image before segmentation; 
        if used, must specify all 6 limits
**Format**: quoted list of 6 integers separated by commas in the order zlo,zhi,ylo,yhi,xlo,xhi
default="0,-1,0,-1,0,-1"	

**--min_slice_padding**
**Description**: Sets the minimum width or height of slices used for segmentation;
        typically needed for XZ or YZ slices to meet DeepCell requirements;
        normally a power of two
**Format**: integer
default="512"

**--chunk_size**
**Description**: Specifies how often to report on slices segmented
**Format**: integer
default=100

### Optional argument for balancing or reducing voxel dimensions

**--downsample_vector**
**Description**: Sets whether block-wise downsampling should be done before segmentation;
         e.g., if X,Y,Z voxel dimensions are 0.5,0.5,1 microns, perhaps use '2,2,1' so
        that voxels are cubic and all directions are equivalent; 
       resulting segmentation masks are upsampled to the original resolution before saving
**Format**: quoted list of 3 integers separated by commas in the order Z, Y, X
default="1,1,1"

### Optional arguments for reducing number of slices to segment

**--sampling_interval**
**Description**: Specifies how frequently to segment slices in XY, XZ and YZ directions
        (segmentations for slices that are skipped are duplicated using neighboring slices;
        goal is to reduce total computation time; sampling interval would normally be higher
        for directions with more slices which would normally be the XZ and YZ directions;
        e.g., if X,Y,Z voxel dimensions are 0.1,0.1,0.4 microns, perhaps use '1,4,4' so that
        slices are segmented less frequently in X)
**Format**: quoted list of 3 integers separated by commas in the order XY, XZ, YZ
default='1,1,1'

**--sampling_reduce*
**Description**: Divisor to reduce sampling_interval by and resentment if quality_threshold is not 
        reached; segmentations for previously segmented slices are kept and reused
**Format**: integer
default=2

**--max_tries**
**Description**: Number of times to reduce the sampling intervals before stopping;
         Stops automatically when sampling interval reaches 1,1,1
**Format**: type=int
default=10

**--quality_threshold**
**Description**: Stops decreasing sampling interval if quality score is greater than this value
**Format**: floating point number
default=infinity

### Optional arguments for reducing compute time

**--JI_range**
**Description**: Change the range of values of the Jaccard index that are tried for optimizing detection of overlap between cells in adjacent 2D slices;
          more intervals mean more time in the "Matching 2D cells" phase
**Format**: quoted list of three floating point numbers (start,end,increment)
default='0.0,0.4,5

**--skip_eval
**Description: Sets whether to call CellSegmentationEvaluator for each JI value segmentation;
         setting to True may be useful to reduce compute time (e.g., when debugging)
**Format**: True or False
default=False

**--skip_blender**
**Description: Sets whether to create Blender files for 3D visualization of final segmentation
         Setting to True may be useful to reduce compute time or if not needed
**Format**: True or False
default=False

**--clear_cache**
**Description: 3DCellComposer saves intermediate files to save time when rerunning or resuming
         an interrupted run.  Setting this option to True deletes those files to start from
        scratch (necessary if significant changes in options are made)
**Format**: True or False
default=False

### Optional arguments for optimizing 2D segmentation

**--maxima_threshold**
**Description**: Sets the corresponding DeepCell parameter (lower values favor more cells)
**Format**: floating point number
default=0.075 (which is the DeepCell default)

**--interior_threshold**
**Description**: Sets the corresponding DeepCell parameter (lower values favor fewer cells)
**Format**: floating point number
default=0.2

**--compartment**
**Description**: Sets the corresponding DeepCell parameter ("both", "whole-cell", "nuclear")
           (Used to ignore one of the input channels, if desired)
**Format**: string
default="both"


### Example Command
Here is an example of a complete command:

```bash
python run_3DCellComposer.py ./data/3D_IMC_image.ome.tiff "Ir191" "In115,Y89,Tb159" "La139,Pr141,Eu151,Gd160,Dy162" --segmentation_method "deepcell" --results_path "myresults"
```

Note that the example image file 3D_IMC_image.ome.tiff is too large to be stored with normal github.  It is therefore stored with git lfs which places a pointer to the actual file into the github repo.  In order to retrieve the file use "git lfs pull" - see <https://git-lfs.com> and
<https://graphite.dev/guides/how-to-use-git-lfs-pull> for more information.

In this command, the script processes the image located at `./data/3D_IMC_image.ome.tiff`, utilizing Ir191 as the nuclear marker, the sum of In115, Y89, Tb159 as the cytoplasmic markers, and the sum of La139, Pr141, Eu151, Gd160, Dy162 as the cell membrane markers. The script employs DeepCell as the 2D segmentation model, and the results are stored in the specified folder within the current folder.

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

'run_3DCellComposerUsingPackage.py' is an example python script that calls ThreeDCellComposer. 

## Contact

Robert F. Murphy - murphy@cmu.edu\
Haoran Chen - haoranchcmu@gmail.com

