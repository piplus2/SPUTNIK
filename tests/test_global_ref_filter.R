# Test global reference filter

library(testthat)
library(SPUTNIK)

test_that("global reference filter", {
  x <- bladderMALDIRompp2010(verbose = TRUE)
  mz <- attr(x, "mass")
  shape <- attr(x, "size")

  msX <- msiDataset(values = x, mz = mz, rsize = shape[1], csize = shape[2])
  msX <- normIntensity(msX, "median")
  msX <- varTransform(msX, "log2")
  
  refImg <- refImageContinuous(msX, method = "sum")
  roiImg <- refImageBinaryOtsu(refImg)

  cat("similarity: pearson
")
  gpfPearson <- globalPeaksFilter(
    msiData = msX, referenceImage = refImg,
    method = "pearson", threshold = 0
  )
  cat("similarity: spearman
")
  gpfSpearman <- globalPeaksFilter(
    msiData = msX, referenceImage = refImg,
    method = "spearman", threshold = 0
  )
  cat("similarity: ssim
")
  gpfSSIM <- globalPeaksFilter(
    msiData = msX, referenceImage = refImg,
    method = "ssim",
    threshold = 0
  )
  cat("similarity: nmi
")
  gpfNMI <- globalPeaksFilter(
    msiData = msX, referenceImage = roiImg,
    method = "nmi", threshold = 0
  )

  expect_is(gpfPearson, "list")
  expect_equal(attr(gpfPearson, "peak.filter"), T)
  expect_equal(attr(gpfPearson, "filter"), "globalPeaks")
  expect_is(gpfSpearman, "list")
  expect_equal(attr(gpfSpearman, "peak.filter"), T)
  expect_equal(attr(gpfSpearman, "filter"), "globalPeaks")
  expect_is(gpfSSIM, "list")
  expect_equal(attr(gpfSSIM, "peak.filter"), T)
  expect_equal(attr(gpfSSIM, "filter"), "globalPeaks")
  expect_is(gpfNMI, "list")
  expect_equal(attr(gpfNMI, "peak.filter"), T)
  expect_equal(attr(gpfNMI, "filter"), "globalPeaks")

  expect_equal(length(gpfPearson$sel.peaks), 502)
  expect_equal(length(gpfSpearman$sel.peaks), 586)
  expect_equal(length(gpfSSIM$sel.peaks), 1125)
  expect_equal(length(gpfNMI$sel.peaks), 515)
})
