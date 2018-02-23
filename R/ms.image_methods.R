## set generics for ms.image-class methods
if (is.null(getGeneric("plot")))
  setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
if (is.null(getGeneric("invertImage")))
  setGeneric("invertImage", function(object, ...) standardGeneric("invertImage"))
if (is.null(getGeneric("smoothImage")))
  setGeneric("smoothImage", function(object, ...) standardGeneric("smoothImage"))
if (is.null(getGeneric("closeImage")))
  setGeneric("closeImage", function(object, ...) standardGeneric("closeImage"))
if (is.null(getGeneric("binOtsu")))
  setGeneric("binOtsu", function(object, ...) standardGeneric("binOtsu"))

#' Visualize a MS image.
#' \code{plot} extends the generic function to \link{ms.image-class} objects.
#'
#' @param x \link{ms.image-class} object. See \link{msImage}.
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
          function(x)
          {
            # Are you plotting the binary mask?
            is.bin <- .isBinary(x)

            df <- melt(x@values)

            gg <- ggplot(df, aes(x = df$X1, y = df$X2,
                                 fill = {if(is.bin)
                                 {factor(df$value)} else {df$value}})) +
              geom_raster() +
              xlab("X") + ylab("Y") +
              {if (is.bin)
              {
                scale_fill_grey(start = 0, end = 1)
              } else
              {
                scale_fill_viridis(option = "inferno")
              }} +
              coord_fixed() +
              {if (length(x@name) != 0)
              {
                ggtitle(x@name)
              }} +
              guides(fill = guide_legend(title = "value")) +
              theme_bw()

            plot(gg)
          }
)

#' Invert the colors of a MS image.
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
setMethod(f = "invertImage",
          signature = signature(object = "ms.image"),
          definition = function(object)
          {
            if (.isBinary(object))
            {
              object@values <- (object@values == 0) * 1
            } else
            {
              object@values <- max(object@values) - object@values
            }
            object
          }
)

#' Apply Gaussian smoothing to a MS image.
#'
#' @param object \link{ms.image-class} object. See \link{msImage}.
#' @param sigma numeric (default = 2). Standard deviation of Gaussian kernel.
#'
#' @return \link{ms.image-class} object after applying smoothing.
#'
#' @importFrom spatstat blur
#' @importFrom spatstat as.im
#'
#' @example R/examples/msImage_smoothImage.R
#'
#' @export
#' @aliases smoothImage
#'
setMethod(f = "smoothImage",
          signature = signature(object = "ms.image"),
          definition = function(object, sigma = 2)
          {
            if (sigma == 0) {
              return(object)
            }
            if (sigma < 0) {
              stop("'sigma' must be positive.")
            }

            object@values <- as.matrix(blur(as.im(object@values), sigma = sigma))
            if (object@scaled)
              object@values <- object@values / max(object@values)

            object
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
#' @importFrom autothresholdr auto_thresh
#' @aliases binOtsu
#'
setMethod(f = "binOtsu",
          signature = signature(object = "ms.image"),
          definition = function(object)
          {
            if (.isBinary(object)) {
              bw <- msImage(object@values, "ROI")
            }
            im.size <- dim(object@values)
            im8bit <- cut(x = c(object@values), breaks = 256, labels = F) - 1
            lev <- auto_thresh(im8bit, method = "Otsu")
            values <- (im8bit > lev) * 1
            bw <- msImage(matrix(values, im.size[1], im.size[2]), "ROI")
            bw
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
setMethod(f = "closeImage",
          signature = signature(object = "ms.image"),
          definition = function(object, kern.size = 5)
          {
            if (!.isBinary(object))
            {
              stop("closeImage can be applied on binary images only.")
            }
            object@values <- as.matrix(mclosing_square(as.cimg(object@values), kern.size))
            object
          }
)
