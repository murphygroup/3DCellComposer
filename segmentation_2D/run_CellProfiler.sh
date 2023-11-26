#!/bin/bash
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

source ~/anaconda3/etc/profile.d/conda.sh
conda activate CellProfiler
bash $script_dir/wrapper/get_CellProfiler_cppipe.sh $1 $2 $3
cellprofiler -r -c -p $1/CellProfiler_config.cppipe -i $1 -o $1
rm $1/CellProfiler_config.cppipe
conda deactivate
