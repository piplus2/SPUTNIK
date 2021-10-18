## .refImage
#' @importFrom stats median prcomp
#' @import irlba
.refImage <- function(msiData,
                      method = "sum",
                      mzQuery = numeric(),
                      mzTolerance = Inf,
                      useFullMZ = TRUE,
                      smoothIm = FALSE,
                      smoothSigma = 2,
                      sampleReference = "detected",
                      invertAlign = FALSE,
                      verbose = TRUE) {
  accept.method.ref <- c("sum", "median", "mean", "pca")
  if (!any(method %in% accept.method.ref)) {
    stop("valid method values are: ", paste0(accept.method.ref, collapse = ", "), ".")
  }

  .stopIfNotValidMSIDataset(msiData)

  use.mz.query <- FALSE
  if (length(mzQuery) != 0 && any(is.finite(mzQuery))) {
    use.mz.query <- TRUE
  }
  
  if (use.mz.query && !useFullMZ) {
    stop("Set either mzQuery of useFullMZ = TRUE.")
  }
  if (use.mz.query && length(mzTolerance) == 0) {
    stop("mzTolerance missing.")
  }

  # Match the peaks indices
  if (use.mz.query) {
    mz.indices <- .mzQueryIndices(mzQuery, msiData@mz, mzTolerance, verbose)
  } else if (useFullMZ) {
    mz.indices <- seq(1, length(msiData@mz))
  }

  # Calculate the reference values
  if (length(mz.indices) == 1) {
    ref.values <- msiData@matrix[, mz.indices]
  } else {
    msiData@matrix[msiData@matrix == 0] <- NA
    ref.values <- switch(
      method,

      "sum" = apply(msiData@matrix[, mz.indices], 1, sum, na.rm = T),

      "median" = {
        msiData@matrix[msiData@matrix == 0] <- NA
        apply(msiData@matrix[, mz.indices], 1, median, na.rm = T)
      },

      "mean" = {
        msiData@matrix[msiData@matrix == 0] <- NA
        apply(msiData@matrix[, mz.indices], 1, mean, na.rm = T)
      },

      "pca" = {
        msiData@matrix[is.na(msiData@matrix)] <- 0
        message("Calculating first principal component...\n")
        pca <- prcomp_irlba(msiData@matrix[, mz.indices], 1, center = TRUE, scale. = TRUE)
        if (cor(pca$x[, 1], apply(msiData@matrix, 1, mean)) < 0) {
          pca$x[, 1] <- (-1) * pca$x[, 1]
        }
        out <- (pca$x[, 1] - min(pca$x[, 1])) / (max(pca$x[, 1]) - min(pca$x[, 1]))
        out
      }
    )
    ref.values[is.na(ref.values)] <- 0
  }

  # Reshape
  ref.values <- matrix(ref.values, msiData@nrow, msiData@ncol)

  # Generate image object
  ref.image <- msImage(ref.values, name = paste0("Ref: ", method), scale = T)
  rm(ref.values)

  if (smoothIm) {
    ref.image <- smoothImage(ref.image, smoothSigma)
  }
  
  # Align to reference image
  if (sampleReference == "detected") {
    if (invertAlign) {
      if (cor(c(ref.image@values), c(msiData@numdetected@values)) > 0)
        ref.image <- invertImage(ref.image)
    } else {
      if (cor(c(ref.image@values), c(msiData@numdetected@values)) < 0)
        ref.image <- invertImage(ref.image)
    }
  } else if (sampleReference == "tic") {
    if (invertAlign) {
      if (cor(c(ref.image@values), c(msiData@totalioncount@values)) > 0)
        ref.image <- invertImage(ref.image)
    } else {
      if (cor(c(ref.image@values), c(msiData@totalioncount@values)) < 0)
        ref.image <- invertImage(ref.image)
    }
  } else {
    stop("Invalid `sampleReference` value.")
  }

  return(ref.image)
}

#' Structural similarity index (SSIM).
#'
#' \code{ssim} returns the value of SSIM between two vectors representing the
#' color intensities of two images.
#'
#' @param x numeric array. Image 1 color intensity array.
#' @param y numeric array. Image 2 color intensity array.
#' @param numBreaks numeric. Number of bins for the color histogram.
#'
#' @return value of SSIM defined between 0 and 1.
#'
#' @details SSIM is an image quality measure, given a reference considered as
#' noise-less image. It can be also used as a perceived similarity measure
#' between images. The images are converted by default in 8bit.
#'
#' @author Paolo Inglese \email{p.inglese14@@imperial.ac.uk}
#'
#' @references Wang, Z., Bovik, A. C., Sheikh, H. R., & Simoncelli, E. P. (2004).
#' Image quality assessment: from error visibility to structural similarity.
#' IEEE transactions on image processing, 13(4), 600-612.
#'
#' @importFrom stats sd cov
#' @export
#'
SSIM <- function(x, y, numBreaks = 256) {
  x <- c(x)
  y <- c(y)

  x <- x / max(x)
  y <- y / max(y)
  x.dig <- cut(as.numeric(x), numBreaks, labels = F) - 1
  y.dig <- cut(as.numeric(y), numBreaks, labels = F) - 1
  rm(x, y)

  C1 <- (0.01 * (numBreaks - 1))^2
  C2 <- (0.03 * (numBreaks - 1))^2

  mux <- mean(x.dig)
  muy <- mean(y.dig)
  sigxy <- cov(x.dig, y.dig)
  sigx <- var(x.dig)
  sigy <- var(y.dig)

  ssim <- ((2 * mux * muy + C1) * (2 * sigxy + C2)) / ((mux**2 + muy**2 + C1) * (sigx + sigy + C2))
  stopifnot(ssim >= -1 && ssim <= 1)

  return(ssim)
}

#' Normalized mutual information (NMI).
#'
#' \code{NMI} returns the normalized mutual information between two \code{ms.image}
#' objects. The normalized mutual information is calculated as the mutual information
#' divided by square-root of the product of the entropies. This function makes
#' use of the functions available in \code{infotheo} R package.
#'
#' @param x numeric array. Image 1 color intensity array.
#' @param y numeric array. Image 2 (binary mask).
#' @param numBins numeric. Number of bins for discretizing the image colors.
#'
#' @return NMI value between 0 and 1.
#'
#' @author Paolo Inglese \email{p.inglese14@@imperial.ac.uk}
#'
#' @references Meyer, P. E. (2009). Infotheo: information-theoretic measures.
#' R package. Version, 1(0).
#'
#' @importFrom infotheo mutinformation entropy
#' @export
#'
NMI <- function(x, y, numBins = 256) {
  x.dig <- cut(as.numeric(c(x)), breaks = numBins, labels = F) - 1
  y.dig <- cut(as.numeric(c(y)), breaks = 2, labels = F) - 1
  rm(x, y)

  mi <- mutinformation(x.dig, y.dig, method = "emp")
  ## Add the sign to the mutual information
  mi <- mi * sign(mean(x.dig[y.dig == 1]) - mean(x.dig[y.dig == 0]))

  h1 <- entropy(x.dig, method = "emp")
  h2 <- entropy(y.dig, method = "emp")

  I <- mi / sqrt(h1 * h2)
  stopifnot(I >= -1 && I <= 1)

  return(I)
}
