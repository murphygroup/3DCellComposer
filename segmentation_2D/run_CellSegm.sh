#!/bin/bash
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
repo_dir=`pwd`
cd $script_dir/CellSegm/examples
echo $script_dir/CellSegm/examples
matlab -nodesktop -r "run_CellSegm $repo_dir/$1 $2 $3; quit"
python $script_dir/wrapper/convert_to_indexed_image.py $repo_dir/$1 CellSegm

