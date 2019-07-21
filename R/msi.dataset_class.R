#' \link{msi.dataset-class} S4 class definition containing the information about
#' the mass spectrometry imaging dataset.
#'
#' @slot matrix the peaks intensity matrix. Rows represent pixels, and columns
#' represent peaks.
#' @slot mz vector of matched m/z values.
#' @slot nrow geometrical shape (number of rows) of image.
#' @slot ncol geometrical shape (number of columns) of image.
#' @slot norm normalization method
#' @slot normoffset numeric offset used for the normalization
#' @slot vartr variance stabilizing transformation
#' @slot vartroffset numeric offset used for the variance stabilizing transformation
#'
#' @name msi.dataset-class
#' @rdname msi.dataset-class
#'
#' @author Paolo Inglese \email{p.inglese14@imperial.ac.uk}
#'
#' @export
#'

setClass(

  "msi.dataset",
  slots = list(
    matrix = "matrix",
    mz = "numeric",
    nrow = "integer",
    ncol = "integer",
    norm = "character",
    vartr = "character",
    normoffset = "numeric",
    vartroffset = "numeric"
  ),

  validity = function(object) {
    if (length(dim(object@matrix)) != 2) {
      return("Intensity matrix must be 2-dimensional.")
    }

    if (any(is.na(object@matrix))) {
      return("Intensity matrix contains NAs.")
    }

    if (any(is.infinite(object@matrix))) {
      return("Intensity matrix contains infinites.")
    }

    if (min(object@matrix) < 0) {
      return("Intensity matrix contains negative values.")
    }

    if (sum(apply(object@matrix, 2, var) == 0) > 0) {
      warning("Some variables are constant.")
    }

    if (length(object@mz) != ncol(object@matrix)) {
      return("M/Z and intensity matrix have incompatible dimensions.")
    }

    return(TRUE)
  }
)
