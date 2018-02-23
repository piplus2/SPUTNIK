## Load package
library("SPUTNIK")

## Create ms.image-class object
msIm <- msImage(values = matrix(rnorm(200), 40, 50), name = "test", scale = TRUE)

## Smooth the image colors
msImSmoothed <- smoothImage(msIm, sigma = 5)
