library(SPUTNIK)

fakeBinImage <- matrix(0, 100, 100)
fakeBinImage[sample(prod(dim(fakeBinImage)), 50)] <- 1

fakeBinMsImage <- msImage(values = fakeBinImage, name = 'ROI', scale = F)

# Remove the objects with a number of connected pixels smaller than 5
fakeBinMsImage <- removeSmallObjects(fakeBinMsImage, threshold = 5)