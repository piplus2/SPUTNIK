#' Constructior for \link{msi.dataset-class} objects.
#'
#' \code{msiDataset} returns a \link{msi.dataset-class} object. It
#' containg information about the matched peaks intensities, the geometrical
#' dimensions of the mass spectral image, the matched m/z values.
#'
#' @param values numeric matrix containing the peaks intensities. Rows represent
#' pixels and columns represent peaks.
#' @param mz array of m/z peaks values.
#' @param rsize geometric shape (number of rows) of image.
#' @param csize geometric shape (number of columns) of image.
#'
#' @return \link{msi.dataset-class} object.
#'
#' @details Function used to construct the main object \code{\link{msi.dataset-class}}.
#' This object contains all the information about peaks intensities (intensity
#' matrix), the geometrical shape of the image (rows, columns), and the vector
#' of the matched m/z values, generated during the peak matching process.
#'
#' @examples
#' ## Load package
#' library("SPUTNIK")
#'
#' ## Create the msi.dataset-class object
#' sz <- c(5, 4)
#' x <- matrix(rnorm(sz[1] * sz[2] * 20), sz[1]*sz[2], 20)
#' mz <- sort(sample(100, ncol(x)))
#' msiX <- msiDataset(x, mz, sz[1], sz[2])
#'
#' @author Paolo Inglese \email{p.inglese14@imperial.ac.uk}
#'
#' @export
#' @importFrom methods new
msiDataset <-  function(values, mz, rsize, csize)
{
  if (ncol(values) != length(mz))
  {
    stop("Incompatible dimensions of m/z vector and intensity matrix.")
  }
  if (nrow(values) != rsize * csize)
  {
    stop("Incompatible rsize and csize values for the provided intensity matrix.")
  }

  object <- new("msi.dataset")

  # Remove attributes and dimnames
  for (n in names(attributes(values))[names(attributes(values)) != "dim"])
    attr(values, n) <- NULL

  object@matrix <- values
  object@mz <- mz
  object@nrow <- as.integer(rsize)
  object@ncol <- as.integer(csize)

  object
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
#' matIm <- matrix(rnorm(200), 40, 50)
#' im <- msImage(values = matIm, name = "random", scale = TRUE)
#'
#' @export
#' @importFrom methods new
#'
msImage <- function(values, name = character(), scale = TRUE)
{
  object <- new("ms.image")

  if (scale) {
    values <- values / max(values)
    object@scaled <- TRUE
  } else {
    object@scaled <- FALSE
  }
  object@values <- values
  object@name <- name

  object
}

#' Genrate a peak filter object.
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
#' values. M/Z values are used to check whether the indices correspond to the
#' matched m/z values in the \code{\link{msi.dataset-class}} object.
#'
#' @author Paolo Inglese \email{p.inglese14@imperial.ac.uk}
#'
#' @examples
#' library("SPUTNIK")
#' mz <- seq(100, 195, 5)
#' mzIdx <- sample(100, 20)
#' names(mzIdx) <- mz
#' peaksFilter <- createPeaksFilter(mzIdx)
#'
#' @seealso \link{applyPeaksFilter-msi.dataset-method}
#'
#' @export
#'
createPeaksFilter <- function(peaksIndices)
{
  if (is.null(names(peaksIndices)))
  {
    warning("Names of peak indices vector elements should match the selected m/z values.")
  }
  l <- list(sel.peaks = peaksIndices)
  attr(l, "peak.filter") <- TRUE
  attr(l, "filter") <- "custom"
  return(l)
}
