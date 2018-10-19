# Test cmplete spatial randomness

test_that("CRS filter", {
  
  x <- bladderMALDIRompp2010(verbose = TRUE)
  mz <- attr(x, 'mass')
  shape <- attr(x, 'size')
  
  msX <- msiDataset(values = x, mz = mz, rsize = shape[1], csize = shape[2])
  refRoi <- refAndROIimages(msX, refMethod = 'sum', roiMethod = 'otsu')
  
  cat('Clark-Evans test...\n')
  csrCE <- CSRPeaksFilter(msX,
                          method = "ClarkEvans",
                          calculateCovariate = FALSE,
                          covMethod = "sum",       # --------------------
                          mzQueryCov = numeric(),  # Covariate arguments
                          mzTolerance = numeric(), #
                          useFullMZCov = TRUE,     #
                          smoothCov = FALSE,       #
                          smoothCovSigma = 2,      #
                          invertCov = FALSE,       #---------------------
                          adjMethod = "bonferroni",
                          returnQvalues = TRUE,
                          plotCovariate = FALSE,
                          verbose = TRUE)
  
  cat('Kolmogorov-Smirnov test...\n')
  csrKS <- CSRPeaksFilter(msX,
                          method = "KS",
                          calculateCovariate = TRUE,
                          covMethod = "sum",       # --------------------
                          mzQueryCov = numeric(),  # Covariate arguments
                          mzTolerance = numeric(), #
                          useFullMZCov = TRUE,     #
                          smoothCov = FALSE,       #
                          smoothCovSigma = 2,      #
                          invertCov = FALSE,       #---------------------
                          adjMethod = "bonferroni",
                          returnQvalues = TRUE,
                          plotCovariate = FALSE,
                          verbose = TRUE)
  
  expect_is(csrCE, 'list')
  expect_equal(attr(csrCE, 'names'), c('p.value', 'q.value'))
  expect_is(csrKS, 'list')
  expect_equal(attr(csrKS, 'names'), c('p.value', 'q.value'))
  
  csrFiltCE <- createPeaksFilter(which(csrCE$q.value < 0.05))
  csrFiltKS <- createPeaksFilter(which(csrKS$q.value < 0.05))
  
  expect_equal(length(csrFiltCE$sel.peaks), 373)
  expect_equal(length(csrFiltKS$sel.peaks), 569)
})
