# 3DCellComposer - A Versatile Pipeline Utilizing 2D Cell Segmentation Methods for 3D Cell Segmentation
Haoran Chen and Robert F. Murphy\
Carnegie Mellon University\
V1.0.0 Nov 17, 2023

We introduce 3DCellComposer, a versatile, open-source software designed as a general solution for 3D cell segmentation. It allows users to choose any existing 2D segmentation model appropriate for their tissue or cell type(s) without requiring any additional training. Moreover, we have enhanced our CellSegmentationEvaluator quality evaluation tool to support 3D images. It allows users to compare and select the most suitable 2D segmentation models for 3D tasks, without the need for human annotations to assess performance.

Reference: Haoran Chen and Robert F. Murphy (2023) 3DCellComposer - A Versatile Pipeline Utilizing 2D Cell Segmentation Methods for 3D Cell Segmentation


## Manual

### Overview
The `3DCellComposer` script processes multiplexed imaging data for cell segmentation. It requires specific inputs including the path to the image file and lists of markers for different cell components (nucleus, cytoplasm, and cell membrane).

### Environment

- Clone this GitHub repository.
- Install Conda, if not already installed.
- Run `conda env create -f environment.yml`. This command will create a new Conda environment 3DCellComposer and install all the necessary packages.


### Command Structure
The script is executed with the following structure:

```bash
python run_3DCellComposer.py [image_path] [nucleus_markers] [cytoplasm_markers] [membrane_markers] [--segmentation_method]
```
### Detailed Input Description

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
   
5. **Segmentation Method (Optional)**:
   - **Description**: Choose the 2D segmentation method.
   - **Format**: String (one of "deepcell", "compare", "custom")
   - **Default**: "deepcell"
   - **Options**:
      - "deepcell" - Uses DeepCell segmentation, which performed the best in our evaluation.
      - "compare" - Compares and selects the best method among 7 different methods.
      - "custom" - Uses a user-provided segmentation method. An empty wrapper is provided.

### Example Command
Here is an example of a complete command using the `3DCellComposer` script:

```bash
python run_3DCellComposer.py /data/3D_IMC_image.ome.tiff "Ir191" "In115,Y89,Tb159" "La139,Pr141,Eu151,Gd160,Dy162" --segmentation_method "deepcell"
```

In this command, the script will process the image located at /data/3D_IMC_image.ome.tiff using Ir191 as the nuclear marker, In115, Y89, Tb159 as the cytoplasmic markers, and La139, Pr141, Eu151, Gd160, Dy162 as the cell membrane markers.

### Notes
- Ensure that there are no spaces between the commas and the marker names.
- The script is designed to handle these inputs and convert them into appropriate lists for processing.


## Contact

Robert F. Murphy - murphy@cmu.edu\
Haoran Chen - hrchen@cmu.edu

