# Test constructors

library(testthat)
library(SPUTNIK)

# msiDataset

imageShape <- c(100, 100)
numPixels <- prod(imageShape)
numPeaks <- 500
mzVector <- sort(sample(seq(150, 1000), numPeaks))
print(is.numeric(mzVector))
print(is.array(mzVector))
randIntensity <- matrix(rnorm(numPixels * numPeaks), numPixels, numPeaks)

test_that("msi.dataset constructor", {
  msX <- msiDataset(
    values = randIntensity, mz = mzVector, rsize = imageShape[1],
    csize = imageShape[2]
  )
  expect_s4_class(msX, "msi.dataset")
})

# msImage

randPixels <- matrix(rnorm(numPixels), imageShape[1], imageShape[2])

test_that("ms.image constructor", {
  msIm <- msImage(values = randPixels, name = "Test", scale = TRUE)
  expect_s4_class(msIm, "ms.image")
})

# peakFilter

numKeep <- 100
keepIdx <- sample(numPeaks, numKeep)
names(keepIdx) <- mzVector[keepIdx]

test_that("peak.filter constructor", {
  filt <- createPeaksFilter(peaksIndices = keepIdx)
  expect_is(filt, "list")
  expect_equal(attr(filt, "peak.filter"), TRUE)
  expect_equal(attr(filt, "filter"), "custom")
})
