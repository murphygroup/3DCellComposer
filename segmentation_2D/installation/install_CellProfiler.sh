#!/bin/bash
conda env create -f $1/CellProfiler.yml
#conda create --name CellProfiler
source ~/anaconda3/etc/profile.d/conda.sh
conda activate CellProfiler_test
pip install cellprofiler==4.0.7
pip install numpy==1.20.0
conda deactivate
cp $1/displaydataonimage.py ~/anaconda3/envs/CellProfiler_test/lib/python3.8/site-packages/cellprofiler/modules
