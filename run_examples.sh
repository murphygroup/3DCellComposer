python run_3DCellComposer.py ./data/cropped_aligned_3D_CODEX_59_ZCYX_T0_T0.ome.tif "DAPI" "VIM,KRT" "CD45,CD34,CD90,PECAM1,CD11B" --segmentation_method "deepcell"

python run_3DCellComposerUsingPackage.py ./data/cropped_aligned_3D_CODEX_59_ZCYX_T0_T0.ome.tif "DAPI" "VIM,KRT" "CD45,CD34,CD90,PECAM1,CD11B" --segmentation_method "deepcell"
