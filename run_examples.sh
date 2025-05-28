# note that the example image file 3D_IMC_image.ome.tiff is too large to be
# stored with normal github.  It is therefore stored with git lfs which places
# a pointer to the actual file into the github repo
# To retrieve the file use "git lfs pull" - see https://git-lfs.com and
# https://graphite.dev/guides/how-to-use-git-lfs-pull for more information

#outputs will be in folder specified by --results_path (IMCdeepcell)
# this takes about six hours on a single cpu with gpu
python run_3DCellComposer.py ./data/3D_IMC_image.ome.tiff "Ir191" "In115,Y89,Tb159" "La139,Pr141,Eu151,Gd160,Dy162" --segmentation_method "deepcell" --chunk_size "10" --skipYZ "True" --min_slice_padding "128" --results_path IMCdeepcell

#outputs will be to the same folder as the input image
# this takes about ten hours on a single cpu with gpu
#uncomment these lines to run
#pip install ThreeDCellComposer
#python run_3DCellComposerUsingPackage.py ./data/3D_IMC_image.ome.tiff "Ir191" "In115,Y89,Tb159" "La139,Pr141,Eu151,Gd160,Dy162" --segmentation_method "deepcell"
