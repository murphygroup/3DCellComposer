#!/bin/bash
rm $1/cellsegm
cd /home/haoranch/projects/HuBMAP/2D-3D/script/cellsegm
cd examples
matlab -nodesktop -r "run_cellsegm $1 $2 $3; quit"
python /home/haoranch/projects/HuBMAP/2D-3D/script/convert_to_indexed_image.py $1 cellsegm

