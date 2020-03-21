## set generics for ms.image-class methods
if (is.null(getGeneric("binOtsu"))) {
  setGeneric("binOtsu", function(object, ...) standardGeneric("binOtsu"))
}

if (is.null(getGeneric("closeImage"))) {
  setGeneric("closeImage", function(object, ...) standardGeneric("closeImage"))
}

if (is.null(getGeneric("invertImage"))) {
  setGeneric("invertImage", function(object, ...) standardGeneric("invertImage"))
}

if (is.null(getGeneric("plot"))) {
  setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
}

if (is.null(getGeneric("removeSmallObjects"))) {
  setGeneric("removeSmallObjects", function(object, ...)
    standardGeneric("removeSmallObjects"))
}

if (is.null(getGeneric("smoothImage"))) {
  setGeneric("smoothImage", function(object, ...) standardGeneric("smoothImage"))
}

#' Visualize an MS image.
#' \code{plot} extends the generic function to \link{ms.image-class} objects.
#'
#' @param x \link{ms.image-class} object. See \link{msImage}.
#' @param palette string. Color palette. See \link{viridis}.
#'
#' @import ggplot2
#' @importFrom viridis scale_fill_viridis
#' @importFrom reshape melt
#'
#' @example R/examples/msImage_plot.R
#'
#' @export
#' @rdname plot
#' @aliases plot
#'
setMethod("plot",
  signature = signature(x = "ms.image", y = "missing"),
  function(x, palette = "inferno") {
    # Are you plotting the binary mask?
    is.bin <- .isBinary(x)

    df <- melt(x@values)

    gg <- ggplot(df, aes(
      x = df$X1, y = df$X2,
      fill = {
        if (is.bin) {
          factor(df$value)
        } else {
          df$value
        }
      }
    )) +
      geom_raster() +
      xlab("X") + ylab("Y") +
      # Use the right palette for continuous or binary valued image
      {
        if (is.bin) {
          scale_fill_grey(start = 0, end = 1)
        } else {
          scale_fill_viridis(option = palette)
        }
      } +
      coord_fixed() + {
        if (length(x@name) != 0) {
          ggtitle(x@name)
        }
      } +
      guides(
        fill = guide_legend(title = "value"),
        plot.title = element_text(hjust = 0.5)
      ) +
      theme_bw()

    plot(gg)
  }
)

#' Invert the colors of an MS image.
#'
#' @param object \link{ms.image-class} object. See \link{msImage}.
#'
#' @return \link{ms.image-class} object after inverting colors.
#'
#' @example R/examples/msImage_invertImage.R
#'
#' @export
#' @aliases invertImage
#'
setMethod(
  f = "invertImage",
  signature = signature(object = "ms.image"),
  definition = function(object) {
    if (.isBinary(object)) {
      object@values <- (object@values == 0) * 1
    } else {
      object@values <- max(object@values) - object@values
    }
    object
  }
)

#' Apply Gaussian smoothing to an MS image.
#'
#' @param object \link{ms.image-class} object. See \link{msImage}.
#' @param sigma numeric (default = 2). Standard deviation of the smoothing
#' Gaussian kernel.
#'
#' @return \link{ms.image-class} smoothed msImage.
#'
#' @importFrom spatstat blur
#' @importFrom spatstat as.im
#'
#' @example R/examples/msImage_smoothImage.R
#'
#' @export
#' @aliases smoothImage
#'
setMethod(
  f = "smoothImage",
  signature = signature(object = "ms.image"),
  definition = function(object, sigma = 2) {
    if (sigma == 0) {
      return(object)
    }
    if (sigma < 0) {
      stop("smoothImage: 'sigma' must be positive.")
    }

    object@values <- as.matrix(blur(as.im(object@values), sigma = sigma))

    if (object@scaled) {
      object@values <- object@values / max(object@values)
    }

    return(object)
  }
)

#' Binarize MS image using Otsu's thresholding.
#'
#' @param object \link{ms.image-class} object. See \link{msImage}.
#'
#' @return \link{ms.image-class} object with binary intensities.
#'
#' @example R/examples/msImage_binOtsu.R
#'
#' @export
#' @import imager
#' @aliases binOtsu
#'
setMethod(
  f = "binOtsu",
  signature = signature(object = "ms.image"),
  definition = function(object) {
    if (.isBinary(object)) {
      bw <- msImage(object@values, "ROI")
    } else {
      km <- kmeans(c(object@values), 2)
      m <- which.min(km$centers)
      thr <- max(c(object@values)[km$cluster == m])
      bw <- (object@values > thr) * 1
      bw <- msImage(as.matrix(bw), "ROI")
    }

    return(bw)
  }
)

#' Apply morphological closing to binary image.
#'
#' @param object \link{ms.image-class} object. See \link{msImage}.
#' @param kern.size numeric. Kernel size.
#'
#' @return \link{ms.image-class} object after closing.
#'
#' @example R/examples/msImage_closeImage.R
#'
#' @export
#' @importFrom imager mclosing_square as.cimg
#' @aliases closeImage
#'
setMethod(
  f = "closeImage",
  signature = signature(object = "ms.image"),
  definition = function(object, kern.size = 5) {
    if (!.isBinary(object)) {
      stop("closeImage can be applied on binary images only.")
    }

    object@values <- as.matrix(mclosing_square(
      im = as.cimg(object@values),
      size = kern.size
    ))

    return(object)
  }
)

#' Remove binary ROI objects smaller than user-defined number of pixels
#'
#' @param object \link{ms.image-class} object. See \link{msImage}.
#' @param threshold numeric. Smallest number of connected pixels.
#' @param border numeric (default = 3). Size of the empty border to add before
#' detecting the connected objects. The border is removed at the end of the
#' process. If `border = 0`, no border is added.
#'
#' @return \link{ms.image-class} object after filtering.
#'
#' @example R/examples/msImage_removeSmallObjects.R
#'
#' @export
#' @import imager
#' @aliases removeSmallObjects
#'
setMethod(
  f = "removeSmallObjects",
  signature = signature(object = "ms.image"),
  definition = function(object, threshold = 5, border = 3) {
    if (!.isBinary(object)) {
      stop("'removeSmallObjects' can be applied on binary images only.")
    }

    if (border < 0) {
      stop("'border' must be positive.")
    }

    roiMat <- object@values == 1

    # Add a border to the ROI image. This can help to identify groups of
    # connected pixels close to the borders.
    if (border > 0) {
      roiMat <- addBorderImage(roiMat, border = border)
    }

    roiMat[roiMat == 0] <- NA

    # Identify the connected components
    roiMat <- as.cimg(roiMat)
    CC <- label(roiMat)
    CC <- as.matrix(CC)

    # Remove the border
    if (border > 0) {
      CC <- remBorderImage(CC, border = border)
    }

    # Filter the connected objects with a number of pixels smaller than
    # threshold
    ux <- unique(c(CC))
    ux <- ux[!is.na(ux)]
    numPixelsObjects <- array(NA, length(ux))
    for (i in 1:length(ux)) {
      numPixelsObjects[i] <- sum(CC == ux[i], na.rm = T)
    }

    ux <- ux[numPixelsObjects >= threshold]

    if (length(ux) == 0) {
      warning("All objects were removed. Returning the original image.")
      return(object)
    }

    # Define the new ROI
    newRoi <- array(NA, prod(dim(object@values)))
    for (i in 1:length(ux))
    {
      newRoi[CC == ux[i]] <- 1
    }
    newRoi <- matrix(newRoi, nrow(object@values), ncol(object@values))
    stopifnot(all(sort(unique(c(newRoi))) == 1))

    newRoi[is.na(newRoi)] <- 0

    object <- msImage(values = newRoi, name = "ROI", scale = F)
    return(object)
  }
)
