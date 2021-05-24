
#############################  CD3  ###############################

#cell surface CD4 receptor image
CD3_cell_surface = data$`26BL`$CD3
#hist(CD_8a_cell_surface)
display(CD3_cell_surface)
CD3_cell_surface.norm = normalize(CD3_cell_surface, inputRange = quantile(CD3_cell_surface, c(0,0.99)))
#hist(CD4_cell_surface.norm)
CD3_cell_surface.norm.blur = gblur(CD3_cell_surface.norm, sigma = 1)
display(CD3_cell_surface.norm.blur)

#cell segentation
cells = rgbImage(green = 1.5*CD3_cell_surface.norm.blur)

#cell surface mask
#get cells that stand out of the moving background window
cmask = thresh(CD3_cell_surface.norm.blur, w=7, h=7, offset=0.05)
#get rid of noises
cmask = opening(cmask, makeBrush(3, shape='disc'))
#apply watershed to further segment cells
cmask = watershed(distmap(cmask), 2)
display(CD3_cell_surface.norm.blur)
display(colorLabels(cmask))

#segment regions with seeds = cmask
ctmask = propagate(CD4_cell_surface.norm.blur, seeds = cmask)

display(colorLabels(ctmask))
#display(paintObjects(cmask, CD4_cell_surface.norm.blur, col = "#ff00ff"))
#paint outline of cell mask on image of cells
segmented = paintObjects(ctmask, cells, col = "#ff00ff")

display(cells)
display(segmented)

#############################  H3  ###############################

#histone3 marker image
histone = data$`26BL`$`Histone-H3`
#hist(CD_8a_cell_surface)
display(histone)
histone.norm = normalize(histone, inputRange = quantile(histone, c(0,0.99)))
#hist(CD4_cell_surface.norm)
histone.norm.blur = gblur(histone.norm, sigma = 1)
display(histone.norm.blur)

#cell segentation
cells = rgbImage(green = 1.5*histone.norm.blur)

#cell mask
#get cells that stand out of the moving background window
cmask = thresh(histone.norm.blur, w=7, h=7, offset=0.05)
#get rid of noises
cmask = opening(cmask, makeBrush(3, shape='disc'))
#apply watershed to further segment cells
cmask = watershed(distmap(cmask), 2)
display(histone.norm.blur)
display(colorLabels(cmask))

#segment regions with seeds = cmask
ctmask = propagate(histone.norm.blur, seeds = cmask)

display(colorLabels(ctmask))
#display(paintObjects(cmask, CD4_cell_surface.norm.blur, col = "#ff00ff"))
#paint outline of cell mask on image of cells
segmented = paintObjects(ctmask, cells, col = "#ff00ff")

display(cells)
display(segmented)

############################   CD 20  #############################
#histone3 marker image
CD20 = data$`26BL`$CD20
#hist(CD_8a_cell_surface)
display(CD20)
CD20.norm = normalize(CD20, inputRange = quantile(CD20, c(0,0.99)))
#hist(CD4_cell_surface.norm)
CD20.norm.blur = gblur(CD20.norm, sigma = 1)
display(CD20.norm.blur)

#cell segentation
cells = rgbImage(green = 1.5*CD20.norm.blur)

#cell mask
#get cells that stand out of the moving background window
cmask = thresh(CD20.norm.blur, w=7, h=7, offset=0.05)
#get rid of noises
cmask = opening(cmask, makeBrush(3, shape='disc'))
#apply watershed to further segment cells
cmask = watershed(distmap(cmask), 2)
display(CD20.norm.blur)
display(colorLabels(cmask))

#segment regions with seeds = cmask
ctmask = propagate(CD20.norm.blur, seeds = cmask)

display(colorLabels(ctmask))
#display(paintObjects(cmask, CD4_cell_surface.norm.blur, col = "#ff00ff"))
#paint outline of cell mask on image of cells
segmented = paintObjects(ctmask, cells, col = "#ff00ff")

display(cells)
display(segmented)

#############################  CD4  ###############################

#cell surface CD4 receptor image
CD3_cell_surface = data$`26BL`$CD3
#hist(CD_8a_cell_surface)
display(CD3_cell_surface)
CD3_cell_surface.norm = normalize(CD3_cell_surface, inputRange = quantile(CD3_cell_surface, c(0,0.99)))
#hist(CD4_cell_surface.norm)
CD3_cell_surface.norm.blur = gblur(CD3_cell_surface.norm, sigma = 1)
display(CD3_cell_surface.norm.blur)

#histone3 marker image
img_histone = data$`26BL`$`Histone-H3` # get marker CD3 from sample 26BL
#hist(img_histone)
display(img_histone)
#display(img_histone, method = "raster")
# trim above 99% and use normalize() to re-scale to [0,1]
img_histone.norm = normalize(img_histone,inputRange = quantile(img_histone,c(0,0.99)))
#hist(img_histone.norm)
display(img_histone.norm)
img_histone.norm.blur = gblur(img_histone.norm, sigma = 1)
display(img_histone.norm.blur)

#cell segentation
cells = rgbImage(green = 1.5*CD3_cell_surface.norm.blur, blue = 1.5*img_histone.norm.blur)

#nuclei mask
nmask = thresh(img_histone.norm.blur, w=5, h=5, offset=0.05)
nmask = opening(nmask, makeBrush(5, shape='disc'))
nmask = watershed(distmap(nmask), 2)
display(img_histone.norm.blur)
display(colorLabels(nmask))

#cell surface(regions) mask
#CD3_cell_surface.norm.blur = gblur(CD3_cell_surface, sigma = 1)
cmask = opening(CD4_cell_surface.norm.blur>0.1, makeBrush(5, shape = 'disc'))
display(cmask)
cmask = propagate(CD4_cell_surface.norm.blur, seeds = nmask, mask = cmask)

display(CD3_cell_surface.norm.blur>0.1)
display(colorLabels(cmask))
#display(paintObjects(cmask, CD4_cell_surface.norm.blur, col = "#ff00ff"))
#paint outline of cell mask on image of cells
segmented = paintObjects(cmask, cells, col = "#ff00ff")

display(cells)
display(segmented)

