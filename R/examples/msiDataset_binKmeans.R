## Load package
library("SPUTNIK")

## Create the msi.dataset-class object
sz <- c(5, 4)
x <- matrix(rnorm(sz[1] * sz[2] * 20), sz[1] * sz[2], 20)
x[x < 0] <- 0
mz <- sort(sample(100, ncol(x)))
msiX <- msiDataset(x, mz, sz[1], sz[2])

## Generate binary mask by applying k-means on the entire dataset
refBin <- binKmeans(msiX)

## Plot the mask
plot(refBin)
