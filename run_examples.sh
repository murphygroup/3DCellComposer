# note that the example image file 3D_IMC_image.ome.tiff is too large to be
# stored with normal github.  It is therefore stored with git lfs which places
# a pointer to the actual file into the github repo
# To retrieve the file use "git lfs pull" - see https://git-lfs.com and
# https://graphite.dev/guides/how-to-use-git-lfs-pull for more information
#
# You can also retrieve 3D IMC images directly from the HubMAP Repository
# using the URLs below (the first URL is the same as the example image in the
# github repository).
#
#https://g-d00e7b.09193a.5898.dn.glob.us/d3130f4a89946cc6b300b115a3120b7a/data/3D_image_stack.ome.tiff
#
#https://g-d00e7b.09193a.5898.dn.glob.us/cd880c54e0095bad5200397588eccf81/data/3D_image_stack.ome.tiff
#
#https://g-d00e7b.09193a.5898.dn.glob.us/a296c763352828159f3adfa495becf3e/data/3D_image_stack.ome.tiff 
#
# These correspond to datasets HBM459.CGSD.533 (883 MB), HBM778.VHHR.349 (570 MB), and HBM387.XZWR.467 (2.05 GB), respectively.
# You can also download datasets by using the HuBMAP Command Line Tool
# (https://docs.hubmapconsortium.org/clt/install-hubmap-clt.html)


# this command processes the full image and takes about eight hours on a single
# cpu with gpu
# outputs will be in folder specified by --results_path (IMCdeepcell)
#python run_3DCellComposer.py ./data/3D_IMC_image.ome.tiff "Ir191" "In115,Y89,Tb159" "La139,Pr141,Eu151,Gd160,Dy162" --segmentation_method "deepcell" --chunk_size "10" --skipYZ "True" --min_slice_padding "128" --results_path IMCdeepcell

# this command uses only part of the image (as specified by the --crop_limits
# option) in order to more rapidly verify proper installation and show the
# results; it takes only a little over an hour
# outputs will be in folder IMCcropdeepcell
python run_3DCellComposer.py ./data/3D_IMC_image.ome.tiff "Ir191" "In115,Y89,Tb159" "La139,Pr141,Eu151,Gd160,Dy162" --segmentation_method "deepcell" --chunk_size "10" --skipYZ "True" --min_slice_padding "128" --results_path IMCcropdeepcell --crop_limits "0,-1,0,256,0,256"

# this command uses the older pip package, and takes around 16 hours.
# outputs will be to a "results" subfolder in the original image folder
# 
#pip install ThreeDCellComposer
#python run_3DCellComposerUsingPackage.py ./data/3D_IMC_image.ome.tiff "Ir191" "In115,Y89,Tb159" "La139,Pr141,Eu151,Gd160,Dy162" --segmentation_method "deepcell"
