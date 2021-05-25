library(EBImage)
if(!require(Rcpp)){install.packages("Rcpp")}
Rcpp::sourceCpp("img_utils.cpp")

file = "..\\CyTOF_data_melanoma.rds"
data = readRDS(file)
sample_list = names(data)

cell_mask <- function(sample, marker){
  # if(!(sample %in% sample_list)){
  #   print(sample, "not in sample list.")
  # }
  # markers = names(data$sample)
  # if(!(marker %in% markers)){
  #   throw(sample, "not in markers of ", sample, " .")
  # }
  #marker image
  img = getElement(getElement(data,sample),marker)
  img.norm = EBImage::normalize(img, inputRange = quantile(img, c(0,0.99)))
  assign("img.norm.blur",gblur(img.norm, sigma = 1),envir=.GlobalEnv)
  
  #cell segentation
  cells = rgbImage(green = 1.5*img.norm.blur)
  
  #cell mask
  #get cells that stand out of the moving background window
  cmask = thresh(img.norm.blur, w=10, h=10, offset=0.05)
  #get rid of noises
  cmask = opening(cmask, makeBrush(3, shape='disc'))
  #apply watershed to further segment cells
  cmask = watershed(distmap(cmask), 2)
  cmask = fillHull(cmask)
  return (cmask)
}
cell_segmentation <- function(sample, marker){
  cmask = cell_mask(sample, marker)
  #segment regions with seeds = cmask
  ctmask = propagate(img.norm.blur, seeds = cmask)
  display(colorLabels(ctmask))
  #paint outline of cell mask on image
  cells = rgbImage(green = 1.5*img.norm.blur)
  segmented = paintObjects(ctmask, cells, col = "#ff00ff")
  
  path = paste(".\\",sample,"\\",sep="")
  name = paste(path,sample,"_",marker,".png",sep="")
  writeImage(segmented, name, quality=85)
  display(segmented)
}
cell_contour <- function(sample, marker){
  cmask = cell_mask(sample, marker)
  Red = as.Image(cpp_obj_contour(cmask)>0)
  Red = fillHull(Red)
  Green = img.norm.blur
  segmented = paintObjects(Red, toRGB(Green), col=c("red", "yellow"), opac=c(1, 0), thick=TRUE)
  
  path = paste(".\\",sample,"\\",sep="")
  name = paste(path,sample,"_",marker,"_contour",".png",sep="")
  writeImage(segmented, name, quality=85)
  display(segmented)
}

#show cells having mean > 300, size > 30
# histone.mean = cpp_obj_mean(mseg = cmask, histone)
# histone.size = cpp_obj_size(mseg = cmask)
# display(histone.mean>300 & histone.size >30)

#excluding the cells having size > 50, mean > 500
# histone.norm = normalize(histone, inputRange = quantile(histone,c(0,0.99)))
# display(histone.norm)
# histone.select = (histone.size>50 & histone.mean > 500)
# display(histone.norm * img.select)
# 

