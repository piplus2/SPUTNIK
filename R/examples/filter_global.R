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

## Generate the reference image and the ROI mask
ref.roi <- refAndROIimages(msiX, refMethod = "sum", roiMethod = "otsu")

## Perform global peaks filter
glob.peaks <- globalPeaksFilter(
  msiData = msiX, referenceImage = ref.roi$Reference,
  method = "pearson", threshold = 0
)

## Apply the filter
msiX <- applyPeaksFilter(msiX, glob.peaks)

## Print the new dimensions
print(dim(getIntensityMat(msiX)))
