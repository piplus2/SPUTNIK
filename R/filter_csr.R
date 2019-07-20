#' Performs the peak selection based on complete spatial randomness test.
#'
#' \code{CSRPeaksFilter} returns the significance for the null hypothesis that the
#' spatial distribution of the peak intensities follow a random pattern. A
#' significant p-value (q-values can be returned after applying multiple testing
#' correction) allows to reject the hypothesis that the spatial distribution of
#' a peak signal is random. The tests are performed using the functions available
#' in the \code{statspat} R package.
#'
#' @param msiData \link{msi.dataset-class} object. See \link{msiDataset}.
#' @param method string (default = \code{"ClarkEvans"}). CSR statistical test
#' applied to the peaks signal. Accepted values are:
#' \itemize{
#'    \item "ClarkEvans": performs a test based on the Clark and Evans aggregation
#'    R index. This test evaluates the compares of the nearest-neighbors distances
#'    to the case of purely random pattern.
#'    \item "KS": performs a test of goodness-of-fit between the signal pixels
#'    associated point process pattern and a spatial covariate using the
#'    Kolmogorov-Smirnov test. The covariate is defined by the reference image.
#' }
#'
#' @param covariateImage \link{ms.image-class} object. An image used as covariate
#' (necessary for Kolmogorov-Smirnov test). If NULL, the covariate image is
#' calculated using the method defined by `covMethod`.
#' @param covMethod string (default = \code{"sum"}). Method used to calculate the
#' reference image. Read only when \code{method = "KS"}. Possible values
#' are described in \code{'refAndROIimages'}.
#' @param mzQueryCov numeric. Values of m/z used to calculate the reference image.
#' 2 values are interpreted as interval, multiple or single values are searched
#' in the m/z vector. It should be left unset when using \code{useFullMZCov = TRUE}.
#' Read only when \code{method = "KS"}.
#' @param mzTolerance numeric. Tolerance in PPM to match the \code{mzQueryCov}
#' values in the m/z vector. It should be left unset when using
#' \code{useFullMZCov = TRUE}.Read only when \code{method = "KS"}.
#' @param useFullMZCov logical (default = \code{TRUE}). Whether all the peaks should be
#' used to calculate the covariate image. Read only when \code{method = "KS"}.
#' @param smoothCov logical (default = \code{FALSE}). Whether the covariate image
#' should be smoothed using a Gaussian kernel. Read only when \code{method = "KS"}.
#' @param smoothCovSigma numeric (default = 2). Standard deviation of the smoothing
#' Gaussian kernel. Read only when \code{method = "KS"}.
#' @param invertCov logical (default = \code{FALSE}). Whether the covariate image
#' colors should be inverted.]
#'
#' @param adjMethod string (default = \code{"bonferroni"}). Multiple testing correction
#' method. Possible values coincide with those of the \code{stats::p.adjust} function.
#' @param returnQvalues logical (default = \code{TRUE}). Whether the computed
#' q-values should be returned together with the p-values.
#' @param plotCovariate logical (default = \code{FALSE}). Whether the covariate image
#' should be visualized. Read only when \code{method = "KS"}.
#' @param verbose logical (defaul = \code{TRUE}). Additional output texts are
#' generated.
#' @param ... additional parameters compatible with the \code{statspat} functions.
#' See \link[spatstat]{cdf.test} for "KS" and \link[spatstat]{clarkevans.test}.
#' for "ClarkEvans"
#'
#' @author Paolo Inglese \email{p.inglese14@imperial.ac.uk}
#'
#' @references Baddeley, A., & Turner, R. (2005). Spatstat: an R package for
#' analyzing spatial point patterns. Journal of statistical software, 12(6), 1-42.
#' @references Clark, P.J. and Evans, F.C. (1954) Distance to nearest neighbour
#' as a measure of spatial relationships in populations. Ecology 35, 445–453.
#' @references Berman, M. (1986) Testing for spatial association between a point
#' process and another stochastic process. Applied Statistics 35, 54–62.
#'
#' @example R/examples/filter_csr.R
#' @export
#' @importFrom stats p.adjust.methods p.adjust cor
#' @importFrom spatstat as.im
#'
CSRPeaksFilter <- function(msiData,
                           method = "ClarkEvans",
                           covariateImage = NULL,
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
                           verbose = TRUE,
                           ...) {
  .stopIfNotValidMSIDataset(msiData)
  if (!is.null(covariateImage)) {
    .stopIfNotValidMSImage(covariateImage)
    if (prod(dim(covariateImage@values)) != nrow(msiData@matrix)) {
      stop("incompatible dimensions between the 'msiData' and the provided
           'covariateImage'")
    }
  }

  # Check the statistical method
  accept.methods <- c("KS", "ClarkEvans")
  if (!any(method %in% accept.methods)) {
    stop(
      "CSRPeaksFilter: valid values for 'method' are: ",
      paste0(accept.methods, collapse = ", "), "."
    )
  }

  if (!any(adjMethod %in% p.adjust.methods)) {
    stop(
      "CSRPeaksFilter: accepted values for 'adjMethod' are: ",
      paste0(p.adjust.methods, collapse = ", "), "."
    )
  }

  # Calculate the reference image. This is used as reference for Kolmogorov-
  # Smirnov test.
  if (method == "KS" && is.null(covariateImage)) {
    if (verbose) {
      cat("Calculating the reference image. This may take a while...\n")
    }
    # Convert the reference image for the spatstat test
    covariateImage <- .refImage(
      msiData = msiData,
      method = covMethod,
      mzQuery = mzQueryCov,
      mzTolerance = mzTolerance,
      useFullMZ = useFullMZCov,
      smoothIm = smoothCov,
      smoothSigma = smoothCovSigma,
      invertIm = invertCov,
      verbose = TRUE
    )
    covariateImage@name <- "Covariate"

    if (plotCovariate) {
      plot(covariateImage)
    }
  }

  if (!is.null(covariateImage)) {
    ref.covariate <- as.im(t(covariateImage@values))
  }

  # Scale in [0, 1]
  if (verbose) {
    cat("Scaling all the peaks intensities in [0, 1]\n")
  }
  msiData <- .scale.all(msiData)

  # Calculate the p-value for the ion images
  cat("Starting CSR tests. This may take a while...\n")

  p_ <- array(NA, length(msiData@mz), dimnames = list(msiData@mz))
  for (ion in 1:length(msiData@mz))
  {
    if (verbose && ion %% 500 == 0) {
      cat(ion, "")
    }

    ## Skip constant images
    if (var(msiData@matrix[, ion]) == 0) {
      next()
    }

    # Transform into a 2D matrix
    im <- msImage(matrix(msiData@matrix[, ion], msiData@nrow, msiData@ncol),
      scale = F
    )

    p_[ion] <- .csr.test.im(
      im = im,
      method = method,
      ref.im = ref.covariate,
      ...
    )
  }
  cat("\n")

  # Multiple testing correction
  q_ <- NULL
  if (returnQvalues) {
    q_ <- p.adjust(p = p_, method = adjMethod)
  }

  out <- list(
    p.value = p_,
    q.value = q_
  )

  return(out)
}

## .csr.test.im
#' @importFrom spatstat owin
#' @importFrom spatstat ppp
#' @importFrom spatstat clarkevans.test
#' @importFrom spatstat cdf.test
.csr.test.im <- function(im,
                         method = "ClarkEvans",
                         ref.im = NULL,
                         win = NULL,
                         ...) {
  .stopIfNotValidMSImage(im)

  if (is.null(win)) {
    win <- owin(xrange = c(1, nrow(im@values)), yrange = c(1, ncol(im@values)))
  }

  # Transform the image into a point pattern process
  im.bw <- binOtsu(im)
  pix <- which(im.bw@values == 1, arr.ind = T)

  p <- switch(
    method,
    "ClarkEvans" = {
      im.ppp <- ppp(
        x = pix[, 1], y = pix[, 2],
        marks = c(im@values[im.bw@values == 1]), window = win
      )
      return(clarkevans.test(X = im.ppp, ...)$p.value)
    },
    "KS" = {
      im.ppp <- ppp(x = pix[, 1], y = pix[, 2], window = win)
      return(cdf.test(
        X = im.ppp, covariate = ref.im,
        test = "ks", ...
      )$p.value)
    }
  )
  return(p)
}
