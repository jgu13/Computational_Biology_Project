install.packages("BiocManager")
BiocManager::install("EBImage")
library("EBImage")
file = "..\\CyTOF_data_melanoma.rds"
data = readRDS(file)

cell_surface = data$`21RD`$CD3 # get marker CD3 from sample 21RD
hist(cell_surface)
display(cell_surface)
#image can be displayed using R’s build-in plotting facilities by calling display with the argument method = "raster" 
display(cell_surface, method = "raster")

# Linearly scale the intensity values of an image to a specified range.
# trim above 99% and use normalize() to re-scale to [0,1]
cell_surface.norm = normalize(cell_surface,inputRange = quantile(cell_surface,c(0,0.99)))
hist(cell_surface.norm)
display(cell_surface.norm)
cell_surface.norm.blur = gblur(cell_surface.norm, sigma = 1)
display(cell_surface.norm.blur)

#Otsu's thresholding: reducing a greyscale image to a binary image
threshold = otsu(cell_surface.norm)
cell_surface_thr = cell_surface.norm > threshold
display(cell_surface_thr)

img_histone = data$`21RD`$`Histone-H3` # get marker CD3 from sample 21RD
hist(img_histone)
display(img_histone)
#image can be displayed using R’s build-in plotting facilities by calling display with the argument method = "raster" 
display(img_histone, method = "raster")

# Linearly scale the intensity values of an image to a specified range.
# trim above 99% and use normalize() to re-scale to [0,1]
img_histone.norm = normalize(img_histone,inputRange = quantile(img_histone,c(0,0.99)))
hist(img_histone.norm)
display(img_histone.norm)
img_histone.norm.blur = gblur(img_histone.norm, sigma = 1)
display(img_histone.norm.blur, method = "raster")

#Otsu's thresholding: reducing a greyscale image to a binary image
threshold = otsu(img_histone.norm)
img_histone_norm.thr = img_histone.norm > threshold
display(img_histone_norm.thr, method = "raster")

#adaptive thresholding
mask = thresh(img_histone, w = 10, h = 10, offset = 0.05)
display(mask)
#segment with bwlabel
mask_table1 = bwlabel(mask)
display(colorLabels(mask_table1))

#distmap: contain the distance of each pixel to its nearest pixel
mask_table2 = watershed(distmap(mask), 2)
display(colorLabels(mask_table2))

#Voronoi tessalation: separate regions 
voronoi_table = propagate(seeds = mask, x = mask, lambda = 500)
display(colorLabels(voronoi_table), "raster")

#cell segmentation
cells = rgbImage(green = 1.5*cell_surface.norm.blur, blue = 1.5*img_histone.norm.blur)
display(cells)
#nuclei mask
n_thr = thresh(img_histone.norm.blur, w=5, h=5, offset=0.05)
nmask = watershed(distmap(n_thr), 2)
display(img_histone.norm.blur)
display(nmask)
###
#nmask = opening(nmask, makeBrush(5, shape='disc'))
#display(colorLabels(nmask))
###
#cell surface(regions) mask
#cell_surface.norm.blur = gblur(cell_surface, sigma = 1)
cmask = opening(cell_surface.norm.blur>0.1, makeBrush(5, shape = 'disc'))
cmask = propagate(cell_surface.norm.blur, seeds = nmask, mask = cmask, lambda = 50)

display(cell_surface.norm.blur>0.1)
#display(cmask)
display(paintObjects(cmask, cell_surface.norm.blur, col = "#ff00ff"))
#paint outline of cell mask on image of cells
segmented = paintObjects(cmask, cells, col = "#ff00ff")

display(cells)
display(segmented)
