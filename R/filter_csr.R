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
#' (required for Kolmogorov-Smirnov test).
#' @param adjMethod string (default = \code{"bonferroni"}). Multiple testing correction
#' method. Possible values coincide with those of the \code{stats::p.adjust} function.
#' @param returnQvalues logical (default = \code{TRUE}). Whether the computed
#' q-values should be returned together with the p-values.
#' @param plotCovariate logical (default = \code{FALSE}). Whether the covariate image
#' should be visualized. Read only when \code{method = "KS"}.
#' @param cores integer (default = 1). Number of CPU cores. Parallel computation if
#' greater than 1.
#' @param verbose logical (default = \code{TRUE}). Additional output texts are
#' generated.
#' @param ... additional parameters compatible with the \code{statspat.core} functions.
#' See \link[spatstat.core]{cdf.test} for "KS" and \link[spatstat.core]{clarkevans.test}.
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
#' @import doSNOW foreach
#' @importFrom stats p.adjust.methods p.adjust cor
#' @importFrom spatstat.geom as.im
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
CSRPeaksFilter <- function(msiData,
                           method = "ClarkEvans",
                           covariateImage = NULL,
                           adjMethod = "bonferroni",
                           returnQvalues = TRUE,
                           plotCovariate = FALSE,
                           cores = 1,
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
    stop('`KS` method requires a reference continuous image in `covariateImage`.')
  }
  
  # Convert the reference image for the spatstat test
  
  if (plotCovariate) {
    covariateImage@name <- "Covariate"
    plot(covariateImage)
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
  cat("Testing ion images for complete spatial randomness...\n")
  
  niter <- length(msiData@mz)
  pb <- txtProgressBar(max = niter, style = 3, width = 80)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  if (cores > 1) {
    
    # Run parallel
    
    cl <- makeCluster(cores)
    doSNOW::registerDoSNOW(cl)
    
    i <- NULL
    p_ <- foreach::foreach(i = 1:niter, .combine = c, .options.snow = opts) %dopar% {
      if (var(msiData@matrix[, i]) == 0) {
        return(NA)
      } else {
        im <- msImage(matrix(msiData@matrix[, i], msiData@nrow, msiData@ncol),
                      scale = FALSE)
        return(
          .csr.test.im(im = im, method = method, ref.im = ref.covariate, ...)
        )
      }
    }
    
    stopCluster(cl)
    
  } else {
    
    p_ <- array(NA, length(msiData@mz))
    for (i in 1:length(msiData@mz)) {
      progress(i)
      if (var(msiData@matrix[, i]) == 0) {
        next()
      }
      
      im <- msImage(matrix(msiData@matrix[, i], msiData@nrow, msiData@ncol),
                    scale = FALSE)
      
      p_[i] <- .csr.test.im(
        im = im,
        method = method,
        ref.im = ref.covariate,
        ...
      )
    }
    
  }
  
  names(p_) <- msiData@mz
  close(pb)

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
#' @importFrom spatstat.geom owin
#' @importFrom spatstat.geom ppp
#' @importFrom spatstat.core clarkevans.test
#' @importFrom spatstat.core cdf.test
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
        X = im.ppp, covariate = ref.im, jitter = TRUE,
        test = "ks", ...
      )$p.value)
    }
  )
  return(p)
}
