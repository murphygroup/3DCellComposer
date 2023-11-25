#!/bin/bash

source ~/anaconda3/etc/profile.d/conda.sh
conda activate ACSS_classic
python aics_classic_wrapper.py $1 $2 $3
conda deactivate


