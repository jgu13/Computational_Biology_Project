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
img = ss$`CD20`
img.norm = normalize(img, inputRange = quantile(img, c(0,0.99)))
display(img.norm)
display(rgbImage(green = img.norm))

img1 = ss$`CD4`
img1.norm = normalize(img1, inputRange = quantile(img1, c(0,0.99)))
display(img1.norm)
display(rgbImage(green = img1.norm))

img2 = ss$`CD3`
img2.norm = normalize(img2, inputRange = quantile(img2, c(0,0.99)))
display(img2.norm)
display(rgbImage(green = img2.norm, red = img.norm))

display(rgbImage(green = img2.norm, red = img.norm))


