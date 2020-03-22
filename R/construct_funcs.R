#' Constructor for \link{msi.dataset-class} objects.
#'
#' \code{msiDataset} returns a \link{msi.dataset-class} object. It
#' contains information about the matched peaks intensities, the geometrical
#' dimensions of the mass spectral image, and the common m/z values.
#'
#' @param values numeric matrix containing the peaks intensities. Rows represent
#' pixels and columns represent peaks.
#' @param mz array of m/z peaks values.
#' @param rsize geometric shape (number of rows) of image.
#' @param csize geometric shape (number of columns) of image.
#' @param verbose boolean (default = TRUE). Additional output.
#'
#' @return \link{msi.dataset-class} object.
#'
#' @details Function used to construct the main object \code{\link{msi.dataset-class}}.
#' This object contains all the information about peaks intensities (intensity
#' matrix), the geometrical shape of the image (rows, columns), and the vector
#' of the common m/z values, generated during the peak matching process.
#'
#' @examples
#' ## Load package
#' library("SPUTNIK")
#'
#' ## Create the msi.dataset-class object
#' sz <- c(5, 4)
#' numIons <- 20
#' x <- matrix(rnorm(prod(sz) * numIons), prod(sz), numIons)
#' mz <- sort(sample(100, numIons))
#' msiX <- msiDataset(x, mz, sz[1], sz[2])
#' @author Paolo Inglese \email{p.inglese14@imperial.ac.uk}
#'
#' @export
#' @importFrom methods new
msiDataset <- function(values, mz, rsize, csize, verbose = TRUE) {

  if (verbose) {
    cat("Creating msiDataset object...\n")
  }
  object <- new("msi.dataset", matrix = values, mz = mz, nrow = rsize,
                ncol = csize, norm = "none", vartr = "none",
                normoffset = 0, vartroffset = 0)
  # Remove
  if (verbose) {
    cat("Detecting constant peaks...\n")
  }
  object <- .remove.const.peaks(object)
  
  return(object)
}

#' Constructor for \link{ms.image-class} objects.
#'
#' @param values numeric matrix representing the pixels intensities. Rows and
#' columns represent the geometrical shape of the image.
#' @param name image name.
#' @param scale logical (default = TRUE). Whether the intensities should be
#' scaled in [0, 1].
#'
#' @return \link{ms.image-class} object.
#'
#' @author Paolo Inglese \email{p.inglese14@imperial.ac.uk}
#'
#' @examples
#' ## Load package
#' library("SPUTNIK")
#'
#' ## MS image
#' imShape <- c(40, 50)
#' matIm <- matrix(rnorm(200), imShape[1], imShape[2])
#' im <- msImage(values = matIm, name = "random", scale = TRUE)
#' @export
#' @importFrom methods new
#'
msImage <- function(values, name = character(), scale = TRUE) {

  if (scale) {
    if (min(values) < 0) {
      warnings("Some pixels have negative intensity. Applying minmax scaling\n")
      values <- (values - min(values)) / (max(values) - min(values))
    } else {
      values <- values / max(values)
    }
    scaled <- TRUE
  } else {
    scaled <- FALSE
  }

  object <- new("ms.image", values = values, name = name, scaled = scaled)
  return(object)
}

#' Generate a peak filter object.
#'
#' \link{createPeaksFilter} returns a \code{peak.filter} object.
#'
#' @param peaksIndices a named array representing the selected peaks. Names correspond
#' to the m/z values.
#'
#' @return \code{peak.filter} object.
#'
#' @details Function to create a custom peak that can be subsequently applied using
#' the function \code{\link{applyPeaksFilter-msi.dataset-method}}. Argument of
#' the function is the index vector of the selected peaks named with their m/z
#' values. The m/z values are used to check whether the indices correspond to the
#' common m/z values in the \code{\link{msi.dataset-class}} object.
#'
#' @author Paolo Inglese \email{p.inglese14@imperial.ac.uk}
#'
#' @examples
#' library("SPUTNIK")
#' mz <- seq(100, 195, 5)
#' mzIdx <- sample(100, 20)
#' names(mzIdx) <- mz
#' peaksFilter <- createPeaksFilter(mzIdx)
#' @seealso \link{applyPeaksFilter-msi.dataset-method}
#'
#' @export
#'
createPeaksFilter <- function(peaksIndices) {
  if (is.null(names(peaksIndices))) {
    stop("Names of 'peakIndices' elements should correspond to the M/Z values.")
  }

  l <- list(sel.peaks = peaksIndices)
  attr(l, "peak.filter") <- TRUE
  attr(l, "filter") <- "custom"  # define this as a custom filter

  return(l)
}
