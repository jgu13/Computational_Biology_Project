disc = makeBrush(31, "disc")
disc = disc / sum(disc)
offset = 0.05
nuc_bg = filter2( nuc, disc )
nuc_th = nuc > nuc_bg + offset
display(nuc_th, all=TRUE)

# The distance map, which contains for each pixel the distance 
# to the nearest background pixel, can be obtained by distmap.
img_histone.ws = watershed(distmap(img_histone_norm_thr), 5)

#Voronoi tessalation
voronoiExamp = propagate(seeds = nmask, x = nmask, lambda = 100)
voronoiPaint = colorLabels (voronoiExamp)
display(voronoiPaint)

#floodfill
rgblogo = toRGB(logo)
points = rbind(c(50, 50), c(100, 50), c(150, 50))
colors = c("red", "green", "blue")
rgblogo = floodFill(rgblogo, points, colors)
display( rgblogo )
#Greyscale
display(floodFill(img, rbind(c(200, 300), c(444, 222)), col=0.2, tolerance=0.2))

#Highlighting objects



