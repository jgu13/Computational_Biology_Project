library(EBImage)
img2 = img[1:30,1:30]
img3 = ifelse(img2 > 0, TRUE, FALSE)

img3 = ifelse(img2 > 0, "TRUE", "FALSE")

img3 = ifelse(img2 > 0, "red", "blue")

img4 <- toRGB(img3)
#color conversion
#system command line for png
#a bunch of command lines for zooming
#system()

ss = data$`21RD`
CD20.img = ss$`CD20`
CD20.img.norm = normalize(CD20.img, inputRange = quantile(CD20.img, c(0,0.99)))
display(CD20.img.norm)
display(rgbImage(green = CD20.img.norm))

# img1 = ss$`CD4`
# CD4.img.norm = normalize(img1, inputRange = quantile(img1, c(0,0.99)))
# display(CD4.img.norm)
# display(rgbImage(green = CD4.img.norm))

img2 = ss$`CD3`
CD3.img.norm = normalize(img2, inputRange = quantile(img2, c(0,0.99)))
display(CD3.img.norm)
display(rgbImage(green = CD3.img.norm))

H3 = ss$`Histone-H3`
H3.img.norm = H3.img.norm = normalize(H3, inputRange = quantile(H3, c(0,0.99)))

combined <- rgbImage(green = CD3.img.norm, red = CD20.img.norm)
sample = "21RD"
marker1 = "CD20"
marker2 = "CD3"
path = paste("..\\RawData\\",sample,"\\",sep="")
name = paste(path,marker1,"+",marker2,".tiff",sep="")
display(combined)
writeImage(combined, name)

assign("21RD_CD20+H3",rgbImage( blue = H3.img.norm, red = CD20.img.norm))
sample = "21RD"
marker1 = "H3"
marker2 = "CD20"
path = paste("..\\RawData\\",sample,"\\",sep="")
name = paste(path,marker1,"+",marker2,".tiff",sep="")
display(`21RD_CD20+H3`)
writeImage(`21RD_CD20+H3`, name)

assign("21RD_CD3+H3",rgbImage( blue = H3.img.norm, red = CD3.img.norm))
sample = "21RD"
marker1 = "H3"
marker2 = "CD3"
path = paste("..\\RawData\\",sample,"\\",sep="")
name = paste(path,marker1,"+",marker2,".tiff",sep="")
display(`21RD_CD3+H3`)
writeImage(`21RD_CD3+H3`, name)

# hist(CD3.img.norm)
temp_path <- 'C:\\Users\\admin\\Documents\\temp.tiff'

file1 = file.choose()
img1 = readImage(file1)
display(img1)

file2 = file.choose()
img2 = readImage(file2)
display(img2)

file3 = file.choose()
img3 = readImage(file3)
display(img3)

#take the inverse of probability
CD20_prob_map <- img1
CD3_prob_map <- img2
#combined <- matrix(nrow = nrow(CD20_prob_map),ncol=ncol(CD20_prob_map))

#50% confidence threshold
thresh_CD20_prob <- ifelse(CD20_prob_map > 0.5, CD20_prob_map, 0)
thresh_CD3_prob <- ifelse(CD3_prob_map > 0.5, CD3_prob_map, 0)
#hist(CD20_prob_map)

#red = label for CD20, blue = label for CD3
# display(CD20_prob_map>0.5)
CD20_cells <- ifelse(thresh_CD20_prob > thresh_CD3_prob, "red", "blue")
display(toRGB(CD20_cells))
CD3_cells <- ifelse(thresh_CD3_prob > thresh_CD20_prob, "red", "blue")
display(toRGB(CD3_cells))

# segment CD20_cells using watershed to get objects
watershed(CD20_cells, )
CD20_Fil <- cpp_obj_size

CD20_dir <- "C:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\H3+CD20\\"
CD20_mask<-readImage(paste(CD20_dir,"H3+CD20Test_ProbabilitiesMask.tiff", sep=""))
CD20_prob_map <- readImage(paste(CD20_dir,"Probability for CD20.tiff", sep=""))
CD20_objs <- readImage(paste(CD20_dir, "H3+CD20Test_ProbabilitiesCell_Obj_Image.tiff", sep="")) 

CD3_dir <- "C:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\H3+CD3\\"
CD3_mask<-readImage(paste(CD3_dir,"H3+CD3Test_ProbabilitiesMask.tiff",sep=""))
CD3_prob_map<-readImage(paste(CD3_dir, "Probability for CD3.tiff", sep=""))
CD3_objs<-readImage(paste(CD3_dir, "H3+CD3Test_ProbabilitiesCell_Obj_Image.tiff", sep=""))

#Use a masked, parallel probability map for each cell type to get a probability mean for each object
#mask probability map
CD20_prob_map <- ifelse(CD20_prob_map > 0.5, CD20_prob_map, 0)
CD20_prob_map <- CD20_prob_map*CD20_mask
CD3_prob_map <- ifelse(CD3_prob_map > 0.5, CD3_prob_map, 0)
CD3_prob_map <- CD3_prob_map*CD3_mask

#get all CD20 objects
#convert labeled image to objects list


#count pixels of every object, while also sum up probability from the probability map


#get all CD3 objects

tiff_file <- file.choose()
img<-EBImage::Image(t(tiff::readTIFF(tiff_file, all = FALSE, info = TRUE, native = FALSE, convert = FALSE, indexed = TRUE, as.is = TRUE)))
img<-readImage(tiff_file)

file <- file.choose()
img <- readImage(file)
display(medianFilter(img, size=3))
display(img)
display(medianFilter(img>0.5, size=3))
display(medianFilter(img>0.5, size=1))
display(opening(img>0.5,makeBrush(size=3, 'diamond')))
writeImage(opening(img>0.5,makeBrush(size=3, 'diamond')), "C:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\H3+CD20\\CellMask.tiff")

file <- file.choose()
img <- readImage(file)
display(opening(img>0.5,makeBrush(size=3, 'diamond')))
writeImage(opening(img>0.5,makeBrush(size=3, 'diamond')), "C:\\Users\\admin\\Documents\\mcgill\\CS_and_Biol\\Comp401\\Cell_seg\\RawData\\21RD\\H3+CD3\\CellMask.tiff")



