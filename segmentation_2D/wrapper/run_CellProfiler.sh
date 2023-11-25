#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate CellProfiler
bash /home/haoranch/projects/HuBMAP/2D-3D/script/get_CellProfiler_cppipe.sh $1 $2 $3
cellprofiler -r -c -p $1/CellProfiler_config.cppipe -i $1 -o $1
conda deactivate
