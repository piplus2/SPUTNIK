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

## Calculate the p-values using the Clark Evans test, then apply Benjamini-
## Hochberg correction.
csr <- CSRPeaksFilter(
  msiData = msiX, method = "ClarkEvans",
  calculateCovariate = FALSE, adjMethod = "BH"
)

## Print selected peaks
print(csr$q.value)

## Create a new filter selecting corrected p-values < 0.001
selIdx <- which(csr$q.value < 0.001)
csrFilter <- createPeaksFilter(selIdx)
