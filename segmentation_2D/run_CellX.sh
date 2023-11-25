#!/bin/bash
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
repo_dir=`pwd`

cd $script_dir/CellX
bash get_xml.sh $repo_dir/$1 $2 $3
matlab -nodesktop -r "run_CellX $repo_dir/$1; quit"
cd $repo_dir/$1
rm seeding_00001.png
rm membrane.tif_control.png
rm final_contour1.png
mv final_mask1.png mask_CellX.png
python $script_dir/wrapper/convert_to_indexed_image.py $repo_dir/$1 CellX
