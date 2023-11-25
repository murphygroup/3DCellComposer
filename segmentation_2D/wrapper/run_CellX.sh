#!/bin/bash
rm $1/*CellX*
cd /home/haoranch/projects/HuBMAP/2D-3D/script/CellX
bash get_xml.sh $1 $2 $3
matlab -nodesktop -r "run_CellX $1; quit"
cd $1
rm seeding_00001.png
rm membrane.tif_control.png
rm final_contour1.png
mv final_mask1.png mask_CellX.png
python /home/haoranch/projects/HuBMAP/2D-3D/script/convert_to_indexed_image.py $1 CellX
