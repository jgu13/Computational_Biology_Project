#reticulate package in R only work with 32-bit python
install_miniconda(path="path/to/install/miniconda") #This line will install a 32-bit miniconda which comes with a 32-bit python
# just install_miniconda() will install a miniconda to the default system folder
library(reticulate)

#set an environment
use_condaenv(condaenv='r-reticulate',required=TRUE)
#another way to set an environment
#use_python(python="path/to/miniconda/envs/r-reticulate/python.exe",required = TRUE) 

#to check if the correct environment is used
#py_config() 

#source a python file
source_python("ConsensusCellMask.py")
np <- import("numpy")
#get the segmentation map
#call function segmentation_merge() from sourced python files
prob_mapA_file <- "path\\to\\21RD\\CD20\\21RD_CD20_Orig_Probabilities.npy"
prob_matA <- np$load(prob_mapA_file)

prob_mapB_file = "path\\to\\21RD\\CD3\\21RD_CD3_Orig_Probabilities.npy"
prob_matB = np.load(prob_mapB_file)

seg_mapA_file <- "path\\to\\21RD\\CD20\\Cell_Object_Image.npy"
seg_matA <- np$load(seg_mapA_file)

seg_mapB_file <- "path\\to\\21RD\\CD3\\Cell_Object_Image.npy"
seg_matB <- np$load(seg_mapB_file)

nuc_bin_file = "path\\to\\21RD\\NucleiMask.npy"
nuc_bin = np$load(nuc_bin_file)

seg_maps <- list(seg_matA,seg_matB)
prob_maps <- list(prob_matA,prob_matB)
py$segmentation_merge(seg_maps, prob_maps, nuc_bin)
#You can also use segmentation_merge(seg_maps, prob_maps, nuc_bin, default_dest="path\\to\\save\\segmentationMap image file"). For example:
#segmentation_merge(seg_maps, prob_maps, nuc_bin, default_dest=".\\segmentationMap.png")




