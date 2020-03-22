# Test pixel count filter

library(testthat)
library(SPUTNIK)

test_that("pixel count filter", {
  MIN_NUM_PIXELS <- 9

  x <- bladderMALDIRompp2010(verbose = TRUE)
  mz <- attr(x, "mass")
  shape <- attr(x, "size")

  msX <- msiDataset(values = x, mz = mz, rsize = shape[1], csize = shape[2])
  msX <- normIntensity(msX, "PQN")
  msX <- varTransform(msX, "log2")
  refRoi <- refAndROIimages(msX, refMethod = "sum", roiMethod = "otsu")
  cpfAggr0 <- countPixelsFilter(
    msiData = msX, roiImage = refRoi$ROI,
    minNumPixels = MIN_NUM_PIXELS, aggressive = 0
  )
  cpfAggr1 <- countPixelsFilter(
    msiData = msX, roiImage = refRoi$ROI,
    minNumPixels = MIN_NUM_PIXELS, aggressive = 1
  )
  cpfAggr2 <- countPixelsFilter(
    msiData = msX, roiImage = refRoi$ROI,
    minNumPixels = MIN_NUM_PIXELS, aggressive = 2
  )
  expect_is(cpfAggr0, "list")
  expect_equal(attr(cpfAggr0, "peak.filter"), T)
  expect_equal(attr(cpfAggr0, "filter"), "countPixels")
  expect_is(cpfAggr1, "list")
  expect_equal(attr(cpfAggr1, "peak.filter"), T)
  expect_equal(attr(cpfAggr1, "filter"), "countPixels")
  expect_is(cpfAggr2, "list")
  expect_equal(attr(cpfAggr2, "peak.filter"), T)
  expect_equal(attr(cpfAggr2, "filter"), "countPixels")

  # Number of selected peaks must be: numAggressive0 >= numAggressive1 >=
  # numAggressive >= 2
  expect_true((length(cpfAggr0$sel.peaks) >= length(cpfAggr1$sel.peaks)) &&
    (length(cpfAggr0$sel.peaks) >= length(cpfAggr2$sel.peaks)) &&
    (length(cpfAggr1$sel.peaks) >= length(cpfAggr2$sel.peaks)))
  expect_equal(length(cpfAggr0$sel.peaks), 174)
  expect_equal(length(cpfAggr1$sel.peaks), 152)
  expect_equal(length(cpfAggr2$sel.peaks), 86)
})
