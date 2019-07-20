library(SPUTNIK)

fakeBinImage <- matrix(0, 100, 100)
fakeBinImage[sample(prod(dim(fakeBinImage)), 2000)] <- 1

fakeBinMsImage <- msImage(values = fakeBinImage, name = "ROI", scale = FALSE)

# Remove the objects with a number of connected pixels smaller than 5
fakeBinMsImage <- removeSmallObjects(fakeBinMsImage, threshold = 5)
