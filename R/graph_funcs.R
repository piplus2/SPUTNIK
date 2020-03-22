#' Compute the reference image and the ROI mask.
#'
#' \code{refAndROIimages} returns the reference image, calculated using the
#' \code{refMethod}, and the ROI binary mask, calculated using \code{roiMethod}.
#' These images represent the basic measures for the filters in SPUTNIK.
#'
#' @param msiData \link[SPUTNIK]{msiDataset} object..
#' @param refMethod string (default = "sum"). Method used to calculate the
#' reference image. Valid values are:
#' \itemize{
#'   \item "sum": peak intensities sum
#'   \item "mean": average peak intensities (without zeros)
#'   \item "median": median peak intensities (without zeros)
#'   \item "pca": first principal component scores.
#' }
#' @param roiMethod string (default = "otsu"). Method used to extract the ROI
#' binary mask. Valid values are:
#' \itemize{
#'   \item "otsu": the reference image is binarized using Otsu's thresholding
#'   \item "kmeans": msiData is partitioned in 2 clusters using k-means
#'   \item "kmeans2": k-means is applied with a user-defined number of clusters
#'   (see Details)
#'   \item "supervised": supervised segmentation based on user-defined areas
#'   corresponding to off-sample and sample regions.
#' }
#' @param mzQueryRef numeric. Values of m/z used to calculate the reference image.
#' 2 values are interpreted as interval, multiple or single values are searched
#' in the m/z vector. It overrides the param \code{useFullMZRef}.
#' @param mzTolerance numeric (default = Inf). Tolerance in PPM to match the
#' \code{mzQueryRef} values in the m/z vector.
#' @param useFullMZRef logical (default = TRUE). Whether all the peaks should be
#' used to calculate the reference image.
#' @param smoothRef logical (default = FALSE). Whether the reference image
#' should be smoothed before binarizing. Only valid for \code{roiMethod = "otsu"}.
#' @param smoothSigma numeric (default = 2). Standard deviation of Gaussian
#' kernel.
#' @param invertRef logical (default = FALSE). Whether the reference image
#' colors should be inverted. This can be necessary when the signal is more
#' intense outside the ROI.
#' @param verbose logical (default = TRUE). Additional output text.
#' @param numClusters numeric (default = 4). Only for 'kmeans2' method. Number
#' of clusters.
#' @param sizeKernel 4-D numeric array or numeric (default = 5). Only for 'kmeans2'.
#' Each element of the 4-D array represents the size of the corners square kernels
#' used to determine the off-tissue clusters. The element order is clockwise:
#' top-left, top-right, bottom-left, bottom-right. If negative, the corresponding
#' corner is skipped. If only a single value is passed, the same kernel size is
#' used for the 4 corners.
#' @param numCores numeric (default = 1). Only for 'kmeans2' method. Number of
#' CPU cores for parallel k-means. It must be smaller than the number of
#' available cores.
#'
#' @details Function to extract the reference image from a \code{\link{msi.dataset-class}}
#' object. Two references images are returned, a continuous-valued and a binary-valued.
#' Multiple methods can be used to extract both the continuous and the binary
#' reference images, which afterwards can be used as argument for the \code{\link{globalPeaksFilter}}
#' filter. When 'kmeans2' is applied, the ROI is obtained by merging the sample-related
#' clusters. The user can set a larger number of cluster than 2 (like in 'kmeans'), in
#' such a way a finer segmentation of the sample-related area can be generated.
#' Currently, the off-sample clusters are identified by looking at the most frequent
#' (statistical mode) labels in the corners of the image.
#'
#' @author Paolo Inglese \email{p.inglese14@imperial.ac.uk}
#'
#' @example R/examples/graph_funcs.R
#'
#' @seealso msiDataset, binOtsu, binKmeans
#' @export
#'
refAndROIimages <- function(msiData,
                            refMethod = "sum",
                            roiMethod = "otsu",
                            mzQueryRef = numeric(),
                            mzTolerance = Inf,
                            useFullMZRef = TRUE,
                            smoothRef = FALSE,
                            smoothSigma = 2,
                            invertRef = FALSE,
                            ## Parameters for kmeans2 ##
                            numClusters = 4, # number of clusters
                            sizeKernel = 5, # number of corners pixels used to
                            # identify the off-sample clusters
                            numCores = 1, # parallel computation for k-means2
                            verbose = TRUE) {
  .stopIfNotValidMSIDataset(msiData)

  accept.method.roi <- c("otsu", "kmeans", "kmeans2", "supervised")
  if (!any(roiMethod %in% accept.method.roi)) {
    stop(
      "valid roiMethod values are: ",
      paste0(accept.method.roi, collapse = ", "), "."
    )
  }

  # Reference image
  ref.image <- .refImage(
    msiData = msiData,
    method = refMethod,
    mzQuery = mzQueryRef,
    mzTolerance = mzTolerance,
    useFullMZ = useFullMZRef,
    smoothIm = smoothRef,
    smoothSigma = smoothSigma,
    invertIm = invertRef
  )
  # ROI
  roi.im <- switch(
    roiMethod,
    "otsu" = binOtsu(ref.image),
    "kmeans" = binKmeans(msiData),
    "kmeans2" = binKmeans2(msiData,
      mzQuery = mzQueryRef,
      mzTolerance = mzTolerance,
      useFullMZ = useFullMZRef,
      numClusters = numClusters,
      numCores = numCores,
      kernelSize = sizeKernel,
      verbose = verbose
    ),
    "supervised" = binSupervised(msiData,
      refImage = ref.image,
      mzQuery = mzQueryRef, # Filter m/z values
      useFullMZ = useFullMZRef, #
      mzTolerance = mzTolerance, #
      method = "svm"
    ) # Currently only 'svm' available
  )

  return(list(Reference = ref.image, ROI = roi.im))
}

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
                      invertIm = FALSE,
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
        pca <- prcomp_irlba(msiData@matrix[, mz.indices], 1)
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
  if (invertIm) {
    ref.image <- invertImage(ref.image)
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
