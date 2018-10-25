#' \link{ms.image-class} definition.
#'
#' @slot values numeric 2-D matrix representing the pixel intensity values.
#' @slot name string. Image name used for plotting.
#' @slot scaled logical. Whether the pixels intensities have been scaled in [0, 1]
#' or not.
#'
#' @method \code{\link{msImage}} default
#' @method binOtsu default
#' @method closeImage default
#' @method invertImage default
#' @method plot default
#' @method removeSmallObjects default
#' @method smoothImage default
#'
#' @author Paolo Inglese \email{p.inglese14@imperial.ac.uk}
#'
#' @export
#'
setClass(

  "ms.image",

  slots = list(
    values = "matrix",
    name = "character",
    scaled = "logical"
  ),

  validity = function(object)
  {
    if (length(dim(object@values)) != 2)
    {
      return("values must be 2-D numeric matrix.")
    }

    if (any(is.na(object@values)))
    {
      return("values contain NA")
    }

    if (any(is.infinite(object@values)))
    {
      return("values contains Inf")
    }

    if (min(object@values) < 0 || max(object@values) > 1)
    {
      return("values not between 0 and 1.")
    }

    if (var(object@values) == 0)
    {
      warning("constant values.")
    }
    return(TRUE)
  }

)
