install.packages("BiocManager")
BiocManager::install("EBImage")
library("EBImage")
file <- file.choose()

data = readRDS(file)

img = data$`21RD`$CD3 # get marker CD3 from sample 21RD

hist(img)

EBImage::display(img)

EBImage::display(img, method = "raster")

img.norm = EBImage::normalize(img,inputRange = c(0,quantile(img,0.99))) # trim above 99% and rescale to [0,1]

hist(img.norm)

EBImage::display(img.norm)

img.norm.blur = EBImage::gblur(img.norm, sigma = 1)

EBImage::display(img.norm.blur)
