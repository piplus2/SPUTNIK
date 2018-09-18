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
#' }
#' @param mzQueryRef numeric. Values of m/z used to calculate the reference image.
#' 2 values are interpreted as interval, multiple or single values are searched
#' in the m/z vector. It should be left unset when using \code{useFullMZRef = TRUE}.
#' @param mzTolerance numeric. Tolerance in PPM to match the \code{mzQueryRef}
#' values in the m/z vector. Only valid when \code{useFullMZ = FALSE}.
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
#'
#' @details Function to extract the reference image from a \code{\link{msi.dataset-class}}
#' object. Two references images are returned, a continuous-valued and a binary-valued.
#' Multiple methods can be used to extract both the continuous and the binary
#' reference images, which afterwards can be used as argument for the \code{\link{globalPeaksFilter}}
#' filter.
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
                            mzTolerance = numeric(),
                            useFullMZRef = TRUE,
                            smoothRef = FALSE,
                            smoothSigma = 2,
                            invertRef = FALSE,
                            verbose = TRUE)
{
  accept.method.roi <- c("otsu", "kmeans")
  if (!any(roiMethod %in% accept.method.roi))
  {
    stop("Valid roiMethod values are: ", paste0(accept.method.roi, collapse = ", "), ".")
  }

  # Ref image
  ref.image <- .refImage(msiData = msiData, method = refMethod, mzQuery = mzQueryRef,
                         mzTolerance = mzTolerance, useFullMZ = useFullMZRef,
                         smoothIm = smoothRef, smoothSigma = smoothSigma,
                         invertIm = invertRef)
  # ROI
  roi.im <- switch(roiMethod,
                   "otsu" = binOtsu(ref.image),
                   "kmeans" = binKmeans(msiData))

  return(list(Reference = ref.image, ROI = roi.im))
}

## .refImage
#' @importFrom stats median prcomp
.refImage <- function(msiData,
                      method = "sum",
                      mzQuery = numeric(),
                      mzTolerance = numeric(),
                      useFullMZ = TRUE,
                      smoothIm = FALSE,
                      smoothSigma = 2,
                      invertIm = FALSE,
                      verbose = TRUE)
{
  accept.method.ref <- c("sum", "median", "mean", "pca")
  if (!any(method %in% accept.method.ref))
  {
    stop("Valid method values are: ", paste0(accept.method.ref, collapse = ", "), ".")
  }

  .stopIfNotValidMSIDataset(msiData)

  if (length(mzQuery) == 0 && !useFullMZ)
  {
    stop("mzQuery and useFullMZ are not compatible.")
  }
  if (length(mzQuery) != 0 && useFullMZ)
  {
    stop("mzQuery and useFullMZ are not compatible.")
  }
  if (length(mzQuery) != 0 && length(mzTolerance) == 0)
  {
    stop("mzTolerance missing.")
  }

  # Match the peaks indices
  if (useFullMZ)
  {
    mz.indices <- seq(1, length(msiData@mz))
  } else
  {
    mz.indices <- .mzQueryIndices(mzQuery, msiData@mz, mzTolerance, verbose)
  }

  # Calculate the reference values
  if (length(mz.indices) == 1)
  {
    ref.values <- msiData@matrix[, mz.indices]
  } else
  {
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
        pca <- prcomp(msiData@matrix[, mz.indices])
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

  if (smoothIm)
  {
    ref.image <- smoothImage(ref.image, smoothSigma)
  }
  if (invertIm)
  {
    ref.image <- invertImage(ref.image)
  }

  ref.image
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
SSIM <- function(x, y, numBreaks = 256)
{
  x <- c(x)
  y <- c(y)
  
  x <- x / max(x)
  y <- y / max(y)
  x.dig <- cut(as.numeric(c(x)), numBreaks, labels = F) - 1
  y.dig <- cut(as.numeric(c(y)), numBreaks, labels = F) - 1
  
  C1 <- (0.01 * (numBreaks - 1)) ^ 2
  C2 <- (0.03 * (numBreaks - 1)) ^ 2
  C3 <- C2 / 2
  
  luminance <- function(x, y)
  {
    ux <- mean(x)
    uy <- mean(y)

    return ((2 * ux * uy + C1) / ((ux * ux) + (uy * uy) + C1))
  }
  
  contrast <- function(x, y) {
    sx2 <- var(x)
    sy2 <- var(y)
    sx <- sqrt(sx2)
    sy <- sqrt(sy2)
    
    return ((2 * sx * sy + C2) / (sx2 + sy2 + C2))
  }
  
  structure <- function(x, y)
  {
    sx <- sqrt(var(x))
    sy <-  sqrt(var(y))
    sxy <- cov(x, y)
    
    return ((sxy + C3) / (sx * sy + C3))   
  }
  
  
  ssim <- luminance(x, y) * contrast(x, y) * structure(x, y)
  
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
#' @param y numeric array. Image 2 binary mask.
#' @param numBins numeric. Number of bins for discretizing the image colors. See
#' \link[infotheo]{discretize}.
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
NMI <- function(x, y, numBins = 256)
{
  x.dig <- cut(as.numeric(c(x)), breaks = numBins, labels = F) - 1
  y.dig <- cut(as.numeric(c(y)), breaks = 2, labels = F) - 1

  mi <- mutinformation(x.dig, y.dig, method = "emp")
  ## Add the sign to the mutual information
  mi <- mi * sign(mean(x.dig[y.dig == 1]) - mean(x.dig[y.dig == 0]))
  
  h1 <- entropy(x.dig, method = "emp")
  h2 <- entropy(y.dig, method = "emp")

  return(mi / sqrt(h1 * h2))
}
