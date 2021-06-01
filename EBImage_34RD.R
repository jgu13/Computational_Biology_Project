
#############################  CD3  ###############################

#cell surface CD4 receptor image
CD3 = data$`34RD`$CD3
#hist(CD_8a_cell_surface)
display(CD3)
CD3.norm = normalize(CD3, inputRange = quantile(CD3, c(0,0.99)))
#hist(CD4_cell_surface.norm)
CD3.norm.blur = gblur(CD3.norm, sigma = 1)
display(CD3.norm.blur)

#cell segentation
cells = rgbImage(green = 1.5*CD3.norm.blur)

#cell surface mask
#get cells that stand out of the moving background window
cmask = thresh(CD3.norm.blur, w=7, h=7, offset=0.05)
#get rid of noises
cmask = opening(cmask, makeBrush(3, shape='disc'))
#apply watershed to further segment cells
cmask = watershed(distmap(cmask), 2)
display(CD3.norm.blur)
display(colorLabels(cmask))

#segment regions with seeds = cmask
ctmask = propagate(CD3.norm.blur, seeds = cmask)

display(colorLabels(ctmask))
#display(paintObjects(cmask, CD4_cell_surface.norm.blur, col = "#ff00ff"))
#paint outline of cell mask on image of cells
segmented = paintObjects(ctmask, cells, col = "#ff00ff")

display(cells)
display(segmented)

#############################  H3  ###############################

#histone3 marker image
histone = data$`34RD`$`Histone-H3`
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
cmask = thresh(histone.norm.blur, w=10, h=10, offset=0.05)
#get rid of noises
cmask = opening(cmask, makeBrush(3, shape='disc'))
#apply watershed to further segment cells
cmask = watershed(distmap(cmask), 2)
cmask = fillHull(cmask)
display(histone.norm.blur)
display(bwlabel(cmask))

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
CD20 = data$`34RD`$CD20
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
cmask = fillHull(cmask)
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
