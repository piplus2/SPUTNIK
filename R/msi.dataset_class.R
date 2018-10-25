#' \link{msi.dataset-class} S4 class definition containing the information about
#' the mass spectrometry imaging dataset.
#'
#' @slot matrix the peaks intensity matrix. Rows represent pixels, and columns
#' represent peaks.
#' @slot mz vector of matched m/z values.
#' @slot nrow geometrical shape (number of rows) of image.
#' @slot ncol geometrical shape (number of columns) of image.
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
    ncol = "integer"
  ),

  validity = function(object)
  {
    if (length(dim(object@matrix)) != 2)
    {
      return("values must be 2-D numeric matrix.")
    }

    if (any(is.na(object@matrix)))
    {
      return("values contain NA")
    }

    if (any(is.infinite(object@matrix)))
    {
      return("values contains Inf")
    }

    if (min(object@matrix) < 0)
    {
      return("negative values.")
    }

    if (sum(apply(object@matrix, 2, var) == 0) > 0)
    {
      warning("values constant.")
    }

    if (length(object@mz) != ncol(object@matrix))
    {
      return("mz and intensities incompatible dimensions.")
    }
    
    return(TRUE)
  }
)
