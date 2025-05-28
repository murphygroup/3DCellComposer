#!/bin/bash
BASEDIR=$(dirname "$0")
#source ~/anaconda3/etc/profile.d/conda.sh
#conda activate Cellpose-3.1.1.1
python $BASEDIR/wrapper/Cellpose-3.1.1.1.py $1 $2
#conda deactivate

