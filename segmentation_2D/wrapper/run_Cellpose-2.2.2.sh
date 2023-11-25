#!/bin/bash
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
python_script="$script_dir/Cellpose-2.2.2_wrapper.py"

source ~/anaconda3/etc/profile.d/conda.sh
conda activate Cellpose-2.2.2
python $python_script $1 XY
python $python_script $1 XZ
python $python_script $1 YZ
conda deactivate



