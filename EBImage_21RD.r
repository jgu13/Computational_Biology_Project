if(!require(EBImage)){
  install.packages("BiocManager")
  BiocManager::install("EBImage")
}
library("EBImage")
file = "..\\CyTOF_data_melanoma.rds"
data = readRDS(file)


########################   CD3  #########################

CD3_cell_surface = data$`21RD`$CD3 # get marker CD3 from sample 21RD
hist(CD3_cell_surface)
display(CD3_cell_surface)
#image can be displayed using R’s build-in plotting facilities by calling display with the argument method = "raster" 
display(CD3_cell_surface, method = "raster")

# Linearly scale the intensity values of an image to a specified range.
# trim above 99% and use normalize() to re-scale to [0,1]
CD3_cell_surface.norm = normalize(CD3_cell_surface,inputRange = quantile(CD3_cell_surface,c(0,0.99)))
hist(CD3_cell_surface.norm)
display(CD3_cell_surface.norm)
CD3_cell_surface.norm.blur = gblur(CD3_cell_surface.norm, sigma = 1)
display(CD3_cell_surface.norm.blur)

#Otsu's thresholding: reducing a greyscale image to a binary image
#threshold = otsu(CD3_cell_surface.norm)
#CD3_cell_surface_thr = CD3_cell_surface.norm > threshold
#display(CD3_cell_surface_thr)

#Otsu's thresholding: reducing a greyscale image to a binary image
#threshold = otsu(img_histone.norm)
#img_histone_norm.thr = img_histone.norm > threshold
#display(img_histone_norm.thr, method = "raster")

#adaptive thresholding
#mask = thresh(img_histone, w = 10, h = 10, offset = 0.05)
#display(mask)
#segment with bwlabel
#mask_table1 = bwlabel(mask)
#display(colorLabels(mask_table1))

#distmap: contain the distance of each pixel to its nearest pixel
#mask_table2 = watershed(distmap(mask), 2)
#display(colorLabels(mask_table2))

#Voronoi tessalation: separate regions 
#voronoi_table = propagate(seeds = mask, x = mask, lambda = 500)
#display(colorLabels(voronoi_table), "raster")

#cell segmentation
cells = rgbImage(green = 1.5*CD3_cell_surface.norm.blur)
display(cells)

#cell mask
#get cells that stand out of a moving background window
cmask = thresh(CD3_cell_surface.norm.blur, w=7, h=7, offset=0.05)
#get rid of noises
cmask = opening(cmask, makeBrush(5, shape='disc'))
#apply watershed to further segment cells
cmask = watershed(distmap(cmask), 2)
display(CD3_cell_surface.norm.blur)
display(colorLabels(cmask))

#regions mask
ctmask = propagate(CD3_cell_surface.norm.blur, seeds = cmask)

#display(CD3_cell_surface.norm.blur>0.1)
#display(colorLabels(cmask))
#display(paintObjects(ctmask, CD3_cell_surface.norm.blur, col = "#ff00ff"))

#paint outline of cell mask on image of cells
segmented = paintObjects(ctmask, cells, col = "#ff00ff")

display(cells)
display(segmented)

#######################   histone H3 #####################
img_histone = data$`21RD`$`Histone-H3` # get marker CD3 from sample 21RD
#hist(img_histone)
#display(img_histone)
#image can be displayed using R’s build-in plotting facilities by calling display with the argument method = "raster" 
#display(img_histone, method = "raster")

# Linearly scale the intensity values of an image to a specified range.
# trim above 99% and use normalize() to re-scale to [0,1]
img_histone.norm = normalize(img_histone,inputRange = quantile(img_histone,c(0,0.99)))
#hist(img_histone.norm)
display(img_histone.norm)
img_histone.norm.blur = gblur(img_histone.norm, sigma = 1)
display(img_histone.norm.blur)

#cell segmentation
cells = rgbImage(green = 1.5*img_histone.norm.blur)
display(cells)

#cell mask
#get cells that stand out of a moving background window
cmask = thresh(img_histone.norm.blur, w=7, h=7, offset=0.05)
#get rid of noises
cmask = opening(cmask, makeBrush(5, shape='disc'))
#apply watershed to further segment cells
cmask = watershed(distmap(cmask), 2)
display(img_histone.norm.blur)
display(colorLabels(cmask))

#regions mask
ctmask = propagate(img_histone.norm.blur, seeds = cmask)

#display(img_histone.norm.blur>0.1)
#display(colorLabels(cmask))
#display(paintObjects(ctmask, img_histone.norm.blur, col = "#ff00ff"))

#paint outline of cell mask on image of cells
segmented = paintObjects(ctmask, cells, col = "#ff00ff")

display(cells)
display(segmented)

########################  CD 20 #########################
#obtain smooth cell image
CD20_cell_surface = data$`21RD`$CD20 # get marker CD3 from sample 21RD
hist(CD20_cell_surface)
display(CD20_cell_surface)

# Linearly scale the intensity values of an image to a specified range.
# trim above 99% and use normalize() to re-scale to [0,1]
CD20_cell_surface.norm = normalize(CD20_cell_surface,inputRange = quantile(CD20_cell_surface,c(0,0.99)))
hist(CD20_cell_surface.norm)
display(CD20_cell_surface.norm)
CD20_cell_surface.norm.blur = gblur(CD20_cell_surface.norm, sigma = 1)
display(CD20_cell_surface.norm.blur)

#cell segmentation
cells = rgbImage(green = 1.5*CD20_cell_surface.norm.blur)
display(cells)

#cell mask
#get cells that stand out of a moving background window
cmask = thresh(CD20_cell_surface.norm.blur, w=5, h=5, offset=0.05)
#get rid of noises
cmask = opening(cmask, makeBrush(3, shape='disc'))
#apply watershed to further segment cells
cmask = watershed(distmap(cmask), 2)
display(CD20_cell_surface.norm.blur)
display(colorLabels(cmask))

#regions mask
ctmask = propagate(CD20_cell_surface.norm.blur, seeds = cmask)

#display(CD20_cell_surface.norm.blur>0.1)
#display(colorLabels(cmask))
#display(paintObjects(ctmask, CD20_cell_surface.norm.blur, col = "#ff00ff"))

#paint outline of cell mask on image of cells
segmented = paintObjects(ctmask, cells, col = "#ff00ff")

display(cells)
display(segmented)

########################  CD 8a #########################

CD8a_cell_surface = data$`21RD`$CD8a
#hist(CD_8a_cell_surface)
#display(CD_8a_cell_surface)
CD8a_cell_surface.norm = normalize(CD8a_cell_surface, inputRange = quantile(CD8a_cell_surface, c(0,0.99)))
#hist(CD_8a_cell_surface.norm)
CD8a_cell_surface.norm.blur = gblur(CD8a_cell_surface.norm, sigma = 1)
display(CD8a_cell_surface.norm.blur)

cells = rgbImage(green = 1.5*CD8a_cell_surface.norm.blur, blue = 1.5*img_histone.norm.blur)

#cell surface(regions) mask
#CD8a_cell_surface.norm.blur = gblur(CD8a_cell_surface, sigma = 1)
cmask = opening(CD8a_cell_surface.norm.blur>0.1, makeBrush(5, shape = 'disc'))
display(cmask)
cmask = propagate(CD8a_cell_surface.norm.blur, seeds = nmask, mask = cmask)

display(CD8a_cell_surface.norm.blur>0.1)
display(colorLabels(cmask))
display(paintObjects(cmask, CD8a_cell_surface.norm.blur, col = "#ff00ff"))
#paint outline of cell mask on image of cells
segmented = paintObjects(cmask, cells, col = "#ff00ff")

display(cells)
display(segmented)

########################   CD 4 #########################

CD4_cell_surface = data$`21RD`$CD4
#hist(CD_8a_cell_surface)
display(CD_8a_cell_surface)
CD4_cell_surface.norm = normalize(CD4_cell_surface, inputRange = quantile(CD4_cell_surface, c(0,0.99)))
#hist(CD4_cell_surface.norm)
CD4_cell_surface.norm.blur = gblur(CD4_cell_surface.norm, sigma = 1)
display(CD4_cell_surface.norm.blur)

cells = rgbImage(green = 1.5*CD4_cell_surface.norm.blur, blue = 1.5*img_histone.norm.blur)

#cell surface(regions) mask
#CD8a_cell_surface.norm.blur = gblur(CD8a_cell_surface, sigma = 1)
cmask = opening(CD4_cell_surface.norm.blur>0.1, makeBrush(3, shape = 'disc'))
display(cmask)
cmask = propagate(CD4_cell_surface.norm.blur, seeds = nmask, mask = cmask)

display(CD4_cell_surface.norm.blur>0.1)
display(colorLabels(cmask))
#display(paintObjects(cmask, CD4_cell_surface.norm.blur, col = "#ff00ff"))
#paint outline of cell mask on image of cells
segmented = paintObjects(cmask, cells, col = "#ff00ff")

display(cells)
display(segmented)

########################   CD 31 #########################

CD4_cell_surface = data$`21RD`$CD4
#hist(CD_8a_cell_surface)
display(CD_8a_cell_surface)
CD4_cell_surface.norm = normalize(CD4_cell_surface, inputRange = quantile(CD4_cell_surface, c(0,0.99)))
#hist(CD4_cell_surface.norm)
CD4_cell_surface.norm.blur = gblur(CD4_cell_surface.norm, sigma = 1)
display(CD4_cell_surface.norm.blur)

cells = rgbImage(green = 1.5*CD4_cell_surface.norm.blur, blue = 1.5*img_histone.norm.blur)

#cell surface(regions) mask
#CD31_cell_surface.norm.blur = gblur(CD31_cell_surface, sigma = 1)
cmask = opening(CD4_cell_surface.norm.blur>0.1, makeBrush(3, shape = 'disc'))
display(cmask)
cmask = propagate(CD4_cell_surface.norm.blur, seeds = nmask, mask = cmask)

display(CD4_cell_surface.norm.blur>0.1)
display(colorLabels(cmask))
#display(paintObjects(cmask, CD4_cell_surface.norm.blur, col = "#ff00ff"))
#paint outline of cell mask on image of cells
segmented = paintObjects(cmask, cells, col = "#ff00ff")

display(cells)
display(segmented)



