library(EBImage)
library(reticulate)
np <- import("numpy")
Rcpp::sourceCpp("ConsensusCellMask.cpp")
source_python("ConsensusCellMask.py")

# seg_mapA_file <- "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD20\\Cell_Object_Image.npy"
# seg_matA <- np$load(seg_mapA_file)
# seg_mapA <- Image(t(seg_matA))
# prob_mapA_file <- "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD20\\21RD_CD20_Orig_Probabilities.tiff"
# #prob_matA <- np$load(prob_mapA_file)
# prob_mapA <- readImage(prob_mapA_file)
# seg_mapB_file <- "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD3\\Cell_Object_Image.npy"
# seg_matB <- np$load(seg_mapB_file)
# seg_mapB <- Image(t(seg_matB))
# prob_mapB_file <- "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD3\\21RD_CD3_Orig_Probabilities.tiff"
# #prob_matB <- np$load(prob_mapB_file)
# prob_mapB <- readImage(prob_mapB_file)
# 
# #seg_mapA <- Image(t(tiff::readTIFF(seg_mapA_file, all = FALSE, info = TRUE, native = FALSE, convert = FALSE, indexed = TRUE, as.is = TRUE)))
# #seg_mapB <- Image(t(tiff::readTIFF(seg_mapB_file, all = FALSE, info = TRUE, native = FALSE, convert = FALSE, indexed = TRUE, as.is = TRUE)))
# 
# merged_mat <- cpp_segmentation_merge(seg_mapA, prob_mapA, seg_mapB, prob_mapB)
# merged_mat_colored <- ifelse(merged_mat==2, "#f94552", merged_mat)
# merged_mat_colored <- ifelse(merged_mat==3, "#bada55", merged_mat)
# merged_mat_colored <- ifelse(merged_mat=="0", "#000000", merged_mat)
# merged_img <- Image(merged_mat_colored)
# display(merged_img)
# 
# seg_mapC_file <- "C:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD31\\Cell_Object_Image.npy"
# seg_matC <- np$load(seg_mapC_file)
# seg_mapC <- Image(t(seg_matC))
# prob_mapC_file <- "C:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD31\\21RD_CD31_Orig_Probabilities.tiff"
# #prob_matA <- np$load(prob_mapA_file)
# prob_mapC <- readImage(prob_mapC_file)
# seg_mapD_file <- "C:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD68\\Cell_Object_Image.npy"
# seg_matD <- np$load(seg_mapD_file)
# seg_mapD <- Image(t(seg_matD))
# prob_mapD_file <- "C:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD68\\21RD_CD68_Orig_Probabilities.tiff"
# #prob_matB <- np$load(prob_mapB_file)
# prob_mapD <- readImage(prob_mapD_file)
# 
# #seg_matC <- Image(t(tiff::readTIFF(seg_matC_file, all = FALSE, info = TRUE, native = FALSE, convert = FALSE, indexed = TRUE, as.is = TRUE)))
# #seg_mapB <- Image(t(tiff::readTIFF(seg_mapB_file, all = FALSE, info = TRUE, native = FALSE, convert = FALSE, indexed = TRUE, as.is = TRUE)))
# 
# merged_mat <- cpp_segmentation_merge(seg_mapC, prob_mapC, seg_mapD, prob_mapD)
# merged_mat_colored <- ifelse(merged_mat==2, "#f97654", merged_mat)
# merged_mat_colored <- ifelse(merged_mat==3, "#d0a7c6", merged_mat)
# merged_mat_colored <- ifelse(merged_mat=="0", "#000000", merged_mat)
# merged_img <- Image(merged_mat_colored)
# display(merged_img)


prob_matA_file <- "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD20\\21RD_CD20_Orig_Probabilities.npy"
prob_mapA <- np$load(prob_matA_file)
prob_matA <- matrix(prob_mapA, nrow=200, ncol=260, byrow=TRUE)

prob_matB_file <- "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD3\\21RD_CD3_Orig_Probabilities.npy"
prob_mapB <- np$load(prob_matB_file)
prob_matB <- matrix(prob_mapB, nrow=200, ncol=260, byrow=TRUE)

prob_matC_file <- "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD31\\21RD_CD31_Orig_Probabilities.npy"
prob_mapC <- np$load(prob_matC_file)
prob_matC <- matrix(prob_mapC, nrow=200, ncol=260, byrow=TRUE)

prob_matD_file <- "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD68\\21RD_CD68_Orig_Probabilities.npy"
prob_mapD <- np$load(prob_matD_file)
prob_matD <- matrix(prob_mapD, nrow=200, ncol=260, byrow=TRUE)

prob_matE_file <- "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD4\\21RD_CD4_Orig_Probabilities.npy"
prob_mapE <- np$load(prob_matE_file)
prob_matE <- matrix(prob_mapE, nrow=200, ncol=260, byrow=TRUE)

seg_mapA_file <- "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD20\\Cell_Object_Image.npy"
seg_matA <- np$load(seg_mapA_file)

seg_mapB_file <- "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD3\\Cell_Object_Image.npy"
seg_matB <- np$load(seg_mapB_file)

seg_mapC_file <- "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD31\\Cell_Object_Image.npy"
seg_matC <- np$load(seg_mapC_file)

seg_mapD_file <- "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD68\\Cell_Object_Image.npy"
seg_matD <- np$load(seg_mapD_file)

seg_mapE_file <- "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\CD4\\Cell_Object_Image.npy"
seg_matE <- np$load(seg_mapE_file)

nuc_bin_file = "E:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\TestWithIllumination\\NucleiMask.npy"
nuc_bin = np$load(nuc_bin_file)
nuc_bin[nuc_bin == TRUE]=1
nuc_bin[nuc_bin == FALSE]=0

seg_maps <- list(seg_matA, seg_matB)
prob_maps <- list(prob_matA, prob_matB)
py$segmentation_merge(seg_maps, prob_maps, nuc_bin)
# marker_mat <- list[[1]]
# prob_mat <- list[[2]]
# outline_mat <- list[[3]]
# 
# marker_mat[marker_mat == 0] <- "#000000"
# marker_mat[marker_mat == 1] <- "#f94552"
# marker_mat[marker_mat == 2] <- "#bada55"
# marker_mat[marker_mat == 3] <- "#f97654"
# marker_mat[marker_mat == 4] <- "#d0a7c6"

# display(Image(t(marker_img)))





