## .isMSImageObject
#' @importFrom methods is
.isMSImageObject <- function(x) {
  return(is(object = x, class2 = "ms.image"))
}

## .isMSIDatasetObject
#' @importFrom methods is
.isMSIDatasetObject <- function(x) {
  return(is(object = x, class2 = "msi.dataset"))
}

## .isBinary
.isBinary <- function(x) {
  .stopIfNotValidMSImage(x)
  if (all(unique(c(x@values)) %in% c(0, 1))) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

## .isPeakFilter
.isPeakFilter <- function(x) {
  if (is.null(attr(x, "peak.filter"))) {
    return(FALSE)
  }
  if (!attr(x, "peak.filter")) {
    return(FALSE)
  }
  return(TRUE)
}

## .stopIfNotValidMSImage
.stopIfNotValidMSImage <- function(x) {
  if (!.isMSImageObject(x)) {
    stop("Not a SPUTNIK::ms.image object.")
  }
  return(invisible(NULL))
}

## .stopIfNotValidMSIDataset
.stopIfNotValidMSIDataset <- function(x) {
  if (!.isMSIDatasetObject(x)) {
    stop("Not a SPUTNIK::msi.dataset object.")
  }
  return(invisible(NULL))
}

## .stopIfNotValidPeakFilter
.stopIfNotValidPeakFilter <- function(x) {
  if (!.isPeakFilter(x)) {
    stop("Not a valid peak.filter result.")
  }
  return(invisible(NULL))
}

## .stopIfNotValidGlobalMethod
.stopIfNotValidGlobalMethod <- function(x) {
  accept.values <- c("pearson", "spearman", "ssim", "nmi")
  if (!any(x %in% accept.values)) {
    stop(
      "Accepted method values are: ",
      paste0(accept.values, collapse = ", "), "."
    )
  }
  return(invisible(NULL))
}

## Validity function for m/z vector
.checkValidMZ <- function(mz) {
  if (!is.numeric(mz)) {
    cat("M/Z must be numeric\n")
    return(FALSE)
  }
  if (length(mz) != length(unique(mz))) {
    cat("M/Z contains duplicated values.\n")
    return(FALSE)
  }
  if (any(diff(mz) <= 0)) {
    cat("M/Z values must be sorted.\n")
    return(FALSE)
  }
  return(TRUE)
}

## Validity function for intensity matrix
.checkValidIntensityMatrix <- function(intensityMatrix, showWarnings = FALSE) {
  if (!is.matrix(intensityMatrix) || !is.numeric(intensityMatrix)) {
    cat("Intensity matrix must be a numeric matrix.\n")
    return(FALSE)
  }
  if (length(dim(intensityMatrix)) != 2) {
    cat("Intensity matrix must be 2-dimensional.\n")
    return(FALSE)
  }
  if (any(!is.finite(intensityMatrix))) {
    cat("Intensity matrix must contain finite values.\n")
    return(FALSE)
  }
  if (any(intensityMatrix < 0)) {
    if (showWarnings) {
      cat("WARNING: Intensity matrix contains negative values.\n")
    }
  }
  if (any(intensityMatrix == 0)) {
    if (showWarnings) {
      cat("WARNING: Zeros present in the intensity matrix. A positive offset 
          may be necessary for normalization or variance stabilizing transformation.\n")
    }
  }
  return(TRUE)
}

## Validity function for image shape
.checkValidImageShape <- function(x) {
  if (!is.finite(x)) {
    cat("Shape must be finite.\n")
    return(FALSE)
  }
  if (!is.numeric(x) || length(x) != 1) {
    cat("Shape must be a numeric value.\n")
    return(FALSE)
  }
  if (x <= 0) {
    cat("Shape must be a positive number.\n")
    return(FALSE)
  }
  return(TRUE)
}