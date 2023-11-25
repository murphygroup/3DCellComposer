#!/bin/bash
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

source ~/anaconda3/etc/profile.d/conda.sh
conda activate ACSS_classic
python $script_dir/wrapper/ACSS_classic_wrapper.py $1 $2 $3
conda deactivate


