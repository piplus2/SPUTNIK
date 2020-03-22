## Load package
library("SPUTNIK")

## Mass spectrometry intensity matrix
X <- matrix(rnorm(16000), 400, 40)
X[X < 0] <- 0

## Print original dimensions
print(dim(X))

## m/z vector
mzVector <- seq(600, 900, by = (900 - 600) / 39)

## Read the image size
imSize <- c(20, 20)

## Construct the ms.dataset object
msiX <- msiDataset(X, mzVector, imSize[1], imSize[2])

## Extract the ROI using k-means
ref.roi <- refAndROIimages(
  msiData = msiX, refMethod = "sum",
  roiMethod = "otsu", useFullMZRef = TRUE
)

## Perform count pixels filtering
count.sel <- countPixelsFilter(
  msiData = msiX, roiImage = ref.roi$ROI,
  minNumPixels = 4, aggressive = 1
)

## Apply the filter
msiX <- applyPeaksFilter(msiX, count.sel)

## Print new dimensions
print(dim(getIntensityMat(msiX)))
