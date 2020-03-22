# Test cmplete spatial randomness

library(testthat)
library(SPUTNIK)

test_that("CSR filter", {
  x <- bladderMALDIRompp2010(verbose = TRUE)
  mz <- attr(x, "mass")
  shape <- attr(x, "size")

  msX <- msiDataset(values = x, mz = mz, rsize = shape[1], csize = shape[2])
  msX <- normIntensity(msX, "PQN")
  msX <- varTransform(msX, "log2")
  refRoi <- refAndROIimages(msX, refMethod = "sum", roiMethod = "otsu")

  cat("Clark-Evans test...
")
  csrCE <- CSRPeaksFilter(msX,
    method = "ClarkEvans",
    covariateImage = NULL, # Not necessary for Clark-Evans
    covMethod = "sum", # --------------------
    mzQueryCov = numeric(), # Covariate arguments
    mzTolerance = numeric(), #
    useFullMZCov = TRUE, #
    smoothCov = FALSE, #
    smoothCovSigma = 2, #
    invertCov = FALSE, #---------------------
    adjMethod = "bonferroni",
    returnQvalues = TRUE,
    plotCovariate = FALSE,
    verbose = TRUE
  )

  cat("Kolmogorov-Smirnov test...
")
  csrKS <- CSRPeaksFilter(msX,
    method = "KS",
    covariateImage = NULL, # Calculate the covariate
    covMethod = "sum", # --------------------
    mzQueryCov = numeric(), # Covariate arguments
    mzTolerance = numeric(), #
    useFullMZCov = TRUE, #
    smoothCov = FALSE, #
    smoothCovSigma = 2, #
    invertCov = FALSE, #---------------------
    adjMethod = "bonferroni",
    returnQvalues = TRUE,
    plotCovariate = FALSE,
    verbose = TRUE
  )

  cat("Passing covariate as argument for Kolmogorov-Smirnov test...
")
  csrKS2 <- CSRPeaksFilter(msX,
    method = "KS",
    covariateImage = refRoi$ROI, # Calculate the covariate
    covMethod = "sum", # --------------------
    mzQueryCov = numeric(), # Covariate arguments
    mzTolerance = numeric(), #
    useFullMZCov = TRUE, #
    smoothCov = FALSE, #
    smoothCovSigma = 2, #
    invertCov = FALSE, #---------------------
    adjMethod = "bonferroni",
    returnQvalues = TRUE,
    plotCovariate = FALSE,
    verbose = TRUE
  )

  expect_is(csrCE, "list")
  expect_equal(attr(csrCE, "names"), c("p.value", "q.value"))
  expect_is(csrKS, "list")
  expect_equal(attr(csrKS, "names"), c("p.value", "q.value"))
  expect_is(csrKS, "list")
  expect_equal(attr(csrKS2, "names"), c("p.value", "q.value"))

  csrFiltCE <- createPeaksFilter(which(csrCE$q.value < 0.05))
  csrFiltKS <- createPeaksFilter(which(csrKS$q.value < 0.05))
  csrFiltKS2 <- createPeaksFilter(which(csrKS2$q.value < 0.05))

  expect_equal(length(csrFiltCE$sel.peaks), 370)
  expect_equal(length(csrFiltKS$sel.peaks), 542)
  expect_equal(length(csrFiltKS$sel.peaks), 542)
})
