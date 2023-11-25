#!/bin/bash
conda env create -f $1/AICS_classic.yml
source ~/anaconda3/etc/profile.d/conda.sh
conda activate ACSS_classic
pip install aicssegmentation==0.4.2
conda deactivate