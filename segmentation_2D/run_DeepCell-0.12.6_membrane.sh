#!/bin/bash
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
python_script="$script_dir/DeepCell-0.12.6_wrapper.py"

source ~/anaconda3/etc/profile.d/conda.sh
conda activate deepcell-0.12.6
python $python_script $1 XY
python $python_script $1 XZ
python $python_script $1 YZ
conda deactivate




