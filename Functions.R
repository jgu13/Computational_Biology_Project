library(EBImage)
if(require(Rcpp)){install.packages("Rcpp")}
library(reticulate)
Sys.setenv(RETICULATE_MINICONDA_PATH="E:/Users/admin/miniconda/condabin/conda.bat")
conda_python(envname='r-reticulate',conda='E:/Users/admin/miniconda/condabin/conda.bat')
use_condaenv(condaenv='r-reticulate',required=TRUE)
np <- import("numpy")
plt <- import("matplotlib.pyplot")
mpimg <- import ("matplotlib.image")
Rcpp::sourceCpp("img_utils.cpp")
Rcpp::sourceCpp("ConsensusCellMask.cpp")

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

#store all marker images with marker colored red, DNA colored blue
all_colored_DNA <- function(){
  for (i in seq_along(data)){
    DNA = data[[i]]$`Histone-H3`
    DNA = normalize(DNA, inputRange = quantile(DNA, c(0,0.99)))
    for(j in seq_along(data[[i]])){
      img = data[[i]][[j]]
      sample = sample_list[i]
      marker = names(data[[i]][j])
      img.norm = normalize(img, inputRange = quantile(img, c(0,0.99)))
      #img.norm.blur= gblur(img.norm, sigma = 1)
      path = paste("..\\RawData\\",sample,"\\Colored\\",sep="")
      name = paste(path,sample,"_",marker,"Colored.tiff",sep="")
      writeImage(rgbImage(red = img.norm, blue = DNA), name)
    }
  }
}

#store all marker images with marker colored red, DNA colored blue
all_colored_DNA_reversed <- function(){
  for (i in seq_along(data)){
    DNA = data[[i]]$`Histone-H3`
    DNA = normalize(DNA, inputRange = quantile(DNA, c(0,0.99)))
    for(j in seq_along(data[[i]])){
      img = data[[i]][[j]]
      sample = sample_list[i]
      marker = names(data[[i]][j])
      img.norm = normalize(img, inputRange = quantile(img, c(0,0.99)))
      #img.norm.blur= gblur(img.norm, sigma = 1)
      path = paste("..\\RawData\\",sample,"\\Colored\\",sep="")
      name = paste(path,sample,"_",marker,"Colored_reversed.tiff",sep="")
      writeImage(rgbImage(red = DNA, blue = img.norm), name)
    }
  }
}

Combine_two_masks <- function(seg_mapA, prob_mapA, seg_mapB, prob_mapB){
  #get the consensus color-labeled map
  colorMat <- cpp_segmentation_merge(seg_mapA, prob_mapA, seg_mapB, prob_map)
  #convert the result to a colored image
  colorMat <- ifelse(1, "red", "green")
  RESImg <- image(coloredMat)
  display(coloredMat)
}

save_and_show <- function(var_name, file_name){
  f_name = paste("./",file_name,".png",sep="")
  writeImage(var_name,f_name)
  img = mpimg$imread(f_name)
  plt$imshow(img)
  plt$show()
}

#(1-CD31)*Histone

CD31 <- data$`21RD`$CD31
CD3 <- normalize(data$`21RD`$CD3, inputRange = quantile(data$`21RD`$CD3, c(0,0.99))) 
CD20 <- normalize(data$`21RD`$CD20, inputRange = quantile(data$`21RD`$CD20, c(0,0.99))) 

Histone <- normalize(data$`21RD`$`Histone-H3`, inputRange = quantile(data$`21RD`$`Histone-H3`, c(0,0.99)))
Histone <- Histone*(1-CD31)

Cropped_cells = Image(CD31[1:260,1:200])
save_and_show(Cropped_cells,"Cropped_Cells")
Cropped_DNA = Image(Histone[1:260,1:200])
save_and_show(Cropped_DNA, "Cropped_DNA")

# f = makeBrush(3,shape='diamond')
# f = f / sum(f)
# offset = 0.03
MultDNA <- gblur(Cropped_DNA, sigma=5)
# MultDNA_bg <- filter2(MultDNA, f)
# MultDNA <- MultDNA > MultDNA_bg + offset
MultDNA <- thresh(MultDNA,w=10,h=10, offset=0.05)
MultDNA <- opening(MultDNA, kern=makeBrush(1,shape=('diamond')))
MultDNA <- watershed(distmap(MultDNA), tolerance = 1, ext=1)
save_and_show(MultDNA)
CellMask <- gblur(Cropped_cells, sigma=1.3488)
f = makeBrush(11,shape='disc')
f = f / sum(f)
CellMask <- filter2(CellMask, f)
CellMask <- opening(CD31,makeBrush(5,shape='disc'))
Cells <- propagate(CD31, seeds=MultDNA, mask=CellMask)
img <- rgbImage(green = Histone, blue = CD31)
Cells <- paintObjects(Cells, img, col='#ff00ff')
Cells <- paintObjects(MultDNA, Cells, col='#ffff00')
save_and_show(Cells)

