library(EBImage)
Rcpp::sourceCpp("img_utils.cpp")
# img.ws is segmentation results from watershed
img.ws <- readImage("..\\26BL_ws.png")
# display(img.ws)
# img0 is H3 signal
img0 = data$`21RD`$`Histone-H3`
img.R = as.Image(cpp_obj_contour(img.ws)>0)
img.G = normalize(img0,inputRange = c(0,300))
img3 = rgbImage(img.R, img.G)
display(img3)

img.mean = cpp_obj_mean(mseg = img.ws, img0)
img.size = cpp_obj_size(mseg = img.ws)
display(img.mean>300 & img.size >30)

img.norm = normalize(img0, inputRange = quantile(img0,c(0,0.99)))
display(img.norm)
img.select = (img.size>50 & img.mean > 500)
display(img.norm * img.select)
