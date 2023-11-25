#!/bin/bash
conda env create -f $1/CellProfiler.yml
source ~/anaconda3/etc/profile.d/conda.sh
conda activate CellProfiler
pip install numpy==1.23.0
conda deactivate
cp $1/displaydataonimage.py ~/anaconda3/envs/CellProfiler/lib/python3.8/site-packages/cellprofiler/modules
cp $1/saveimages.py ~/anaconda3/envs/CellProfiler/lib/python3.8/site-packages/cellprofiler/modules
