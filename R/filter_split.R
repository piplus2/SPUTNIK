#' Test for the presence of split peaks.
#'
#' \link{splitPeaksFilter} returns a list of estimated split peak indices. Each
#' element of the list contains an array of the original peak indices that can
#' be merged. The name of the list element is the new m/z value associated with
#' the merged peaks.
#'
#' @param msiData \link{msi.dataset-class} object. See \link{msiDataset}.
#' @param mzTolerance numeric (default = 5). Maximum distance in PPM between the
#' m/z values of two peaks to consider them for merging. See 'Details' section.
#' @param sharedPixelsRatio numeric (default = 0). Maximum fraction of common pixels
#' where the signal of two peaks is different from zero to consider them for
#' merging. See 'Details' section.
#' @param sparseness string (default = \code{"scatter.ratio"}). Method used to estimate
#' the 'scatteredness' of the peak image. See 'Details' section.
#' @param threshold numeric (default = 0.5). Threshold for scatteredness measure
#' to consider peaks for merging. At least one of the merging peaks should have
#' a measure associated with presence of structure.
#' @param returnDetails logical (default = \code{TRUE}). Add details on merged peaks
#' in the results.
#' @param verbose logical (default = \code{TRUE}). Additional output text.
#'
#' @details \link{splitPeaksFilter} determines whether close peaks represent the
#' same signal. This estimation is based on multiple conditions:
#'
#'  \enumerate{
#'
#'   \item peaks m/z values should be closer than \code{mzTolerance} PPM
#'
#'   \item at least one of the peak images should be structured, accordingly to
#'   the \code{sparseness} measure. The \code{threshold} determines whether the
#'   pixel images are structured or not. The possible measures are:
#'
#'   \itemize{
#'     \item \code{"scatter.ratio"}: ratio between the number of non-zero pixels
#'     and the image size after binarization using Otsu's thresholding. A value close
#'     to 0 is associated with a more structured image, whereas a value close to
#'     1 is associated with a less structured image. A suggested parameter of
#'     \code{threshold = 0.5} represents the maximum value for this measure for
#'     a structured image. Minimum possible value is 1 / ( # non-zero pixels ).
#'
#'     \item \code{"spatial.chaos"}: similar to the scatter ratio taking into
#'     account of the color histogram. A value close to 1 represents a structured
#'     image, whereas a value close to 0 represents a more scattered image.
#'     A suggested parameter of \code{threshold = 0.8} represents the minimum value
#'     for this measure for a structured image. Maximum possible value is
#'     1 - 1 / ( # histogram bins ). Here, we use the default number of bins
#'     equal to 30.
#'
#'     \item \code{"gini.index"}: Gini index measures the image sparsity. A value
#'     close to 1 is associated with a sparse image whereas a value close to 0
#'     is associated with a more uniform image. A suggested value of \code{threshold = 0.9}
#'     represents the maximum value of this measure for a structured image.
#'   }
#'
#'   \item the merged peaks image should be more structured than the single
#'   peak images, accordingly to the selected \code{sparseness}.
#' }
#'
#' @author Paolo Inglese \email{p.inglese14@imperial.ac.uk}
#'
#' @references Palmer, A., Phapale, P., Chernyavsky, I., Lavigne, R., Fay, D.,
#' Tarasov, A., ... & Becker, M. (2017). FDR-controlled metabolite annotation for
#' high-resolution imaging mass spectrometry. Nature methods, 14(1), 57.
#' @references Hurley, N., & Rickard, S. (2009). Comparing measures of sparsity.
#' IEEE Transactions on Information Theory, 55(10), 4723-4741.
#'
#' @example R/examples/filter_split.R
#'
#' @export
#' @importFrom stats var
#'
splitPeaksFilter <- function(msiData,
                             mzTolerance = 5,
                             sharedPixelsRatio = 0,
                             sparseness = "scatter.ratio",
                             threshold = 0.5,
                             returnDetails = TRUE,
                             verbose = TRUE) {
  .stopIfNotValidMSIDataset(msiData)

  sparse.accept <- c("gini.index", "scatter.ratio", "spatial.chaos")
  if (!any(sparseness %in% sparse.accept)) {
    stop(
      "accepted values for 'sparseness' are: ",
      paste0(sparse.accept, collapse = ", "), "."
    )
  }

  # Determine the peaks closest than mzTolerance
  if (verbose) {
    cat("Determining close peaks. This may take a while...\n")
  }

  grouped.peaks_ <- list()
  k <- 1
  grouped.peaks_[[k]] <- 1
  d <- array(NA, length(msiData@mz))
  for (i in 2:length(msiData@mz))
  {
    # Skip the constant signal
    if (var(msiData@matrix[, i]) == 0) {
      next()
    }
    # Check whether the group of peaks are close enough
    test.peaks_ <- c(grouped.peaks_[[k]], i)
    d[i] <- .distance.ppm.interval(msiData@mz[test.peaks_])
    if (.distance.ppm.interval(msiData@mz[test.peaks_]) <= mzTolerance) {
      # Test whether the peaks signals are not localized in the same spatial
      # regions.
      grouped.bin <- apply(msiData@matrix[, test.peaks_], 2, function(z) z != 0)
      # The sum of the binary images should not contain numbers larger than 1
      grouped.inters <- apply(grouped.bin, 1, sum)
      if (sum(grouped.inters > 1) / length(grouped.inters) <= sharedPixelsRatio) {
        grouped.peaks_[[k]] <- test.peaks_
      } else {
        k <- k + 1
        grouped.peaks_[[k]] <- i
      }
    } else {
      k <- k + 1
      grouped.peaks_[[k]] <- i
    }
  }
  rm(i)

  # Extract the groups of 2 or more peaks
  sel.groups <- which(unlist(lapply(grouped.peaks_, length)) > 1)
  if (verbose) {
    cat("Min. distance between peaks (ppm) =", min(d, na.rm = T), "\n")
  }
  if (length(sel.groups) == 0) {
    cat("No peak groups to merge found.\n")
    return(NULL)
  } else {
    cat("Found", length(sel.groups), "groups of peaks to test for scatteredness.\n")
  }
  grouped.peaks_ <- grouped.peaks_[sel.groups]

  # If the union of the close peaks returns a less scattered image, then merge
  # the peaks into one measure which m/z corresponds to the mean value of the
  # merged m/z values.
  if (verbose) {
    cat("Testing scatteredness of merged images...\n")
  }

  fin.peaks_ <- list()
  k <- 0
  merge.mz <- c()
  for (j in 1:length(grouped.peaks_))
  {
    # Calculate the scatteredness of the ion images and that of the union
    group.peaks_ <- unique(grouped.peaks_[[j]])
    sc_meas <- array(NA, length(group.peaks_) + 1)
    for (l in 1:length(group.peaks_))
    {
      im_ <- matrix(msiData@matrix[, group.peaks_[l]], msiData@nrow, msiData@ncol)
      im_ <- im_ / max(im_)
      sc_meas[l] <- switch(sparseness,
        "scatter.ratio" = scatter.ratio(im = im_),
        "gini.index" = gini.index(x = im_),
        "spatial.chaos" = spatial.chaos(im = im_)
      )
      rm(im_)
    }
    rm(l)

    # Depending on the selected measure, different thresholding are applied
    # Scatter ratio: a structured image has a low value
    # Gini index: a less sparse image has a low value
    # Spatial chaos: a structured image has a high value
    meas.condition <- switch(sparseness,
      "scatter.ratio" = (any(sc_meas[1:(length(sc_meas) - 1)] <= threshold)),
      "gini.index" = (any(sc_meas[1:(length(sc_meas) - 1)] <= threshold)),
      "spatial.chaos" = (any(sc_meas[1:(length(sc_meas) - 1)] >= threshold))
    )
    if (!meas.condition) {
      next()
    }

    im_ <- matrix(apply(msiData@matrix[, group.peaks_], 1, sum), msiData@nrow, msiData@ncol)
    im_ <- im_ / max(im_)
    sc_meas[length(sc_meas)] <- switch(sparseness,
      "scatter.ratio" = scatter.ratio(im = im_),
      "gini.index" = gini.index(x = im_),
      "spatial.chaos" = spatial.chaos(im = im_)
    )
    # Check the regularity after merging
    meas.condition <- switch(sparseness,
      "scatter.ratio" = all(sc_meas[1:(length(sc_meas) - 1)] >= sc_meas[length(sc_meas)]),
      "gini.index" = all(sc_meas[1:(length(sc_meas) - 1)] >= sc_meas[length(sc_meas)]),
      "spatial.chaos" = (any(sc_meas[1:length(sc_meas)] <= threshold))
    )
    if (meas.condition) {
      k <- k + 1
      fin.peaks_[[k]] <- group.peaks_
      merge.mz[k] <- mean(msiData@mz[group.peaks_])
    }
    names(fin.peaks_) <- merge.mz
  }
  if (length(fin.peaks_) == 0) {
    if (verbose) {
      cat("No peak groups to merge found.\n")
    }
    return(NULL)
  }
  if (verbose) {
    cat(length(unlist(fin.peaks_)), "peaks merged in", length(fin.peaks_), "peaks.\n")
  }

  # Return also the list of the m/z values that were merged
  list.merged.mz <- NA
  if (returnDetails) {
    list.merged.mz <- vector(mode = "list", length = length(fin.peaks_))
    for (i in 1:length(fin.peaks_))
    {
      list.merged.mz[[i]] <- msiData@mz[fin.peaks_[[i]]]
    }
    names(list.merged.mz) <- names(fin.peaks_)
  }

  res <- list()
  res$merged.peaks <- fin.peaks_
  res$merged.mz <- list.merged.mz
  res$details <- list.merged.mz
  attr(res, "peak.filter") <- TRUE
  attr(res, "filter") <- "splitPeaks"

  return(res)
}
