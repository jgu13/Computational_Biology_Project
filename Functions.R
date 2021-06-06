library(EBImage)
if(!require(Rcpp)){install.packages("Rcpp")}
Rcpp::sourceCpp("img_utils.cpp")

file = "..\\CyTOF_data_melanoma.rds"
data = readRDS(file)
sample_list = names(data)

cell_mask <- function(sample, marker){
  #marker image
  img = getElement(getElement(data,sample),marker)
  img.norm = normalize(img, inputRange = quantile(img, c(0,0.99)))
  assign("img.norm.blur",gblur(img.norm, sigma = 1),envir=.GlobalEnv)
  path = paste(".\\",sample,"\\",sep="")
  name = paste(path,sample,"_",marker,".png",sep="")
  writeImage(img.norm.blur, name, quality=100)
  
  #cell segentation
  cells = rgbImage(green = 1.5*img.norm.blur)
  
  #cell mask
  #get cells that stand out of the moving background window
  cmask = thresh(img.norm.blur, w=15, h=15, offset=0.05)
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
  writeImage(segmented, name, quality=100)
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
  writeImage(segmented, name, quality=100)
  display(segmented)
}

analyze_all = function(){
  l1<-rep(sample_list,each=9)
  l2<-rep(names(getElement(data,sample_list[1])),times=3)
  mapply(cell_segmentation, l1, l2)
  mapply(cell_contour, l1, l2)
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

#Batch produce normalized cell image
all_data_images <- function(){
  for (i in seq_along(data)){
    for(j in seq_along(data[[i]])){
      img = data[[i]][[j]]
      sample = sample_list[i]
      marker = names(data[[i]][j])
      img.norm = normalize(img, inputRange = quantile(img, c(0,0.99)))
      #img.norm.blur= gblur(img.norm, sigma = 1)
      path = paste("..\\RawData\\",sample,"\\",sep="")
      name = paste(path,sample,"_",marker,".tiff",sep="")
      writeImage(img.norm, name)
    }
  }
}

# display_alt <- function(var_name){
#   
# }



