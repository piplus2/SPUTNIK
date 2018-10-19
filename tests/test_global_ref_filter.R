# Test global reference filter

test_that("global reference filter", {

  x <- bladderMALDIRompp2010(verbose = TRUE)
  mz <- attr(x, 'mass')
  shape <- attr(x, 'size')
  
  msX <- msiDataset(values = x, mz = mz, rsize = shape[1], csize = shape[2])
  refRoi <- refAndROIimages(msX, refMethod = 'sum', roiMethod = 'otsu')
  
  cat('similarity: pearson\n')
  gpfPearson <- globalPeaksFilter(msiData = msX, referenceImage = refRoi$Reference,
                                  method = 'pearson', threshold = 0)
  cat('similarity: spearman\n')
  gpfSpearman <- globalPeaksFilter(msiData = msX, referenceImage = refRoi$Reference,
                                   method = 'spearman', threshold = 0)
  cat('similarity: ssim\n')
  gpfSSIM <- globalPeaksFilter(msiData = msX, referenceImage = refRoi$Reference,
                               method = 'ssim')
  cat('similarity: nmi\n')
  gpfNMI <- globalPeaksFilter(msiData = msX, referenceImage = refRoi$ROI,
                              method = 'nmi', threshold = 0)
  
  expect_is(gpfPearson, 'list')
  expect_equal(attr(gpfPearson, 'peak.filter'), T)
  expect_equal(attr(gpfPearson, 'filter'), 'globalPeaks')
  expect_is(gpfSpearman, 'list')
  expect_equal(attr(gpfSpearman, 'peak.filter'), T)
  expect_equal(attr(gpfSpearman, 'filter'), 'globalPeaks')
  expect_is(gpfSSIM, 'list')
  expect_equal(attr(gpfSSIM, 'peak.filter'), T)
  expect_equal(attr(gpfSSIM, 'filter'), 'globalPeaks')
  expect_is(gpfNMI, 'list')
  expect_equal(attr(gpfNMI, 'peak.filter'), T)
  expect_equal(attr(gpfNMI, 'filter'), 'globalPeaks')
  
  expect_equal(length(gpfPearson$sel.peaks), 560)
  expect_equal(length(gpfSpearman$sel.peaks), 535)
  expect_equal(length(gpfSSIM$sel.peaks), 789)
  expect_equal(length(gpfNMI$sel.peaks), 558)
})
