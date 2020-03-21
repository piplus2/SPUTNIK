#' Filter based on the minimum number of connected pixels in the ROI.
#'
#' \code{countPixelsFilter} selects peaks which signals are localized in regions
#' consisting of a minimum number of connected pixels in the ROI.
#'
#' @param msiData \link{msi.dataset-class} object. See \link{msiDataset}.
#' @param roiImage \link{ms.image-class} object representing the ROI mask. See
#' \link{msImage}.
#' @param minNumPixels integer (default = 9). Smallest number of connected pixels
#' used to select a peak.
#' @param smoothPeakImage logical (default = \code{FALSE}). Whether the peak
#' images should be smoothed before determining the connected components.
#' @param smoothSigma numeric (default = 2). Standard deviation of the smoothing Gaussian
#' kernel.
#' @param closePeakImage logical (default = \code{FALSE}). Whether morphological
#' closing should be applied to the binary peak images.
#' @param closeKernSize numeric (default = 5). Kernel size for the morphological
#' closing operation. Kernel shape is fixed to \code{diamond}.
#' @param aggressive integer (default = 0). Defines the level of aggressiveness
#' of the filter. See 'Details' section.
#' @param verbose logical (default = \code{TRUE}). Additional output text.
#'
#' @details Count filter tries to determine and remove peaks which signal is
#' scattered in a region unrelated with the expected ROI. A minimum number of
#' connected pixels in the ROI is used to trigger the filter. This value should
#' be carefully set equal to the geometrical size of the smallest expected
#' informative sub-region. Each peak image is binarized using Otsu's thresholding
#' and the connected components are extracted. The filter selects those peaks
#' that show, within the ROI, at least one connected component of size larger or
#' equal to \code{minNumPixels}. The level of aggressiveness, associated with
#' increasing values of the parameter \code{aggressive}, determines whether the
#' size of the connected components within the ROI should be compared with that
#' of the connected components localized outside the ROI. If \code{aggressive = 0},
#' no comparison is performed. If \code{aggressive = 1}, the filter checks whether
#' the max size of the connected components localized outside the ROI is smaller
#' or equal to the maximum size of the connected components inside the ROI.
#' If \code{aggressive = 2}, a stricter filter checks whether the maximum size
#' of the connected components localized outside the ROI is smaller than
#' \code{minNumPixels}. Different aggressiveness levels can produce completely
#' different results, depending on the nature of the analyzed dataset.
#'
#' @return \code{peak.filter} object. See \link{applyPeaksFilter-msi.dataset-method}.
#'
#' @author Paolo Inglese \email{p.inglese14@imperial.ac.uk}
#'
#' @example R/examples/filter_count.R
#'
#' @seealso applyPeaksFilter
#' @export
#' @import imager
#'
countPixelsFilter <- function(msiData,
                              roiImage,
                              minNumPixels = 9,
                              smoothPeakImage = FALSE,
                              smoothSigma = 2,
                              closePeakImage = FALSE,
                              closeKernSize = 5,
                              aggressive = 0,
                              verbose = TRUE) {
  .stopIfNotValidMSIDataset(msiData)
  .stopIfNotValidMSImage(roiImage)
  
  accept.aggr.values <- c(0, 1, 2)
  if (!aggressive %in% accept.aggr.values) {
    stop(paste0("aggressive can be equal to ", paste0(accept.aggr.values, collapse = ", "), ".\n"))
  }
  if (minNumPixels >= prod(getShapeMSI(msiData))) {
    stop("minNumPixels must be smaller than number of MSI pixels.")
  }
  if (smoothSigma >= getShapeMSI(msiData)[1] || smoothSigma >= getShapeMSI(msiData)[2]) {
    stop("smoothSigma must be smaller than MSI shape.")
  }
  if (closeKernSize >= min(getShapeMSI(msiData))) {
    stop("closeKernSize must be smaller than MSI shape.")
  }
  if (!is.logical(verbose)) {
    stop("verbose must be logical value.")
  }

  # Count the number of connected pixels within and outside the ROI. In order
  # to accept a peak as informative, there must be at least one group of connected
  # pixels of minNumPixels size within the ROI. This is necessary to discriminate
  # between peaks that are randomly distributed within the ROI.
  if (verbose) {
    cat("Counting connected pixels within signal region...\n")
  }

  filter.results <- array(NA, length(msiData@mz), dimnames = list(msiData@mz))
  max.count.inside <- filter.results
  max.count.outside <- filter.results

  for (i in 1:length(msiData@mz))
  {
    im <- matrix(msiData@matrix[, i], msiData@nrow, msiData@ncol)
    im <- msImage(values = im, name = as.character(msiData@mz[i]), scale = T)
    # Apply smoothing
    if (smoothPeakImage) {
      im <- smoothImage(im, smoothSigma)
    }
    # Binarize using Otsu's thresholding
    im.bw <- binOtsu(im)
    # Morphological closing
    if (closePeakImage) {
      im.bw <- closeImage(im.bw, kern.size = closeKernSize)
    }
    # Count the connected regions inside and outside the ROI
    conn.comps.within <- label(as.cimg(im.bw@values * (roiImage@values == 1)))
    conn.comps.within <- as.matrix(conn.comps.within)
    conn.comps.outside <- label(as.cimg(im.bw@values * (roiImage@values != 1)))
    conn.comps.outside <- as.matrix(conn.comps.outside)

    large.conn <- FALSE

    # Don't consider the background (label = 0)
    conn.table.inside <- table(c(conn.comps.within[conn.comps.within != 0]))
    conn.table.outside <- table(c(conn.comps.outside[conn.comps.outside != 0]))

    if (length(conn.table.inside) == 0) {
      conn.table.inside <- 0
    }
    if (length(conn.table.outside) == 0) {
      conn.table.outside <- 0
    }

    if (aggressive == 0) {
      outside.cond <- TRUE
    } else if (aggressive == 1) {
      # If aggressive = 1, check whether the largest connected component outside
      # the ROI is smaller than the largest connected component within the ROI.
      # This is based on the idea that signal associated with the sample should
      # show more structured patterns inside the ROI.
      outside.cond <- max(conn.table.outside) <= max(conn.table.inside)
    } else if (aggressive == 2) {
      # If aggressive = 2, check if the largest connected component outside the
      # is smaller than minNumPixels. In this way, we are stricter about the possible
      # structuredness of the signal outside the ROI.
      outside.cond <- max(conn.table.outside) < minNumPixels
    }

    if (any(conn.table.inside >= minNumPixels) && outside.cond) {
      large.conn <- TRUE
    }

    max.count.inside[i] <- max(conn.table.inside)
    max.count.outside[i] <- max(conn.table.outside)
    filter.results[i] <- large.conn
  }

  out <- list(
    max.count.inside = max.count.inside,
    max.count.outside = max.count.outside,
    sel.peaks = which(filter.results)
  )
  attr(out, "peak.filter") <- TRUE
  attr(out, "filter") <- "countPixels"

  return(out)
}
