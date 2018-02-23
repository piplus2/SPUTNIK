#' Pixel scatteredness ratio.
#'
#' \code{scatter.ratio} returns a measure of image scatteredness represented by
#' the ratio between the number of connected components and the total number of
#' non-zero pixels. The number of connected components is calculated from the
#' binarized image using Otsu's method.
#'
#' @param im 2-D numeric matrix representing the image pixel intensities.
#'
#' @references Otsu, N. (1979). A threshold selection method from gray-level
#' histograms. IEEE transactions on systems, man, and cybernetics, 9(1), 62-66.
#'
#' @author Paolo Inglese \email{p.inglese14@imperial.ac.uk}
#'
#' @seealso \link{gini.index} \link{spatial.chaos}
#'
#' @example R/examples/sparseness.R
#'
#' @export
#' @importFrom SDMTools ConnCompLabel
#' @importFrom autothresholdr auto_thresh
#'
scatter.ratio <- function(im)
{
  im.size <- dim(im)
  im8bit <- cut(x = c(im), breaks = 256, labels = F) - 1
  lev_ <- auto_thresh(im8bit, method = "Otsu")
  bin_ <- (im8bit > lev_) * 1
  lbl_ <- unique(c(ConnCompLabel(matrix(bin_, im.size[1], im.size[2]))))
  return(length(lbl_[lbl_ != 0]) / sum(bin_))
}

#' Gini index.
#'
#' \code{gini.index} returns the Gini index of the ion intensity vector as a
#' measure of its sparseness. The intensity vector is first quantized in N
#' levels (default = 256). A value close to 1 represents a high level of
#' sparseness, a value close to 0 represents a low level of sparseness.
#'
#' @param x numeric. Peak intensity array.
#' @param levels numeric (default = 256). Number of levels for the peak
#' intensity quantization.
#'
#' @return A value between 0 and 1. High levels of signal sparsity are associated
#' with values close to 1, whereas low levels of signal sparsity are associated
#' with values close to 0.
#'
#' @references Hurley, N., & Rickard, S. (2009). Comparing measures of sparsity.
#' IEEE Transactions on Information Theory, 55(10), 4723-4741.
#'
#' @author Paolo Inglese \email{p.inglese14@imperial.ac.uk}
#'
#' @example R/examples/sparseness.R
#'
#' @seealso \link{scatter.ratio} \link{spatial.chaos}
#' @export
#'
gini.index <- function(x, levels = 256)
{
  x <- cut(x = x, breaks = levels, labels = F) - 1
  h <- sort(x)
  N <- length(h)
  num.zeros <- sum(h == 0)
  h <- h[h != 0]
  l1.norm <- sum(abs(h))
  tmp.term <- 0
  for (i in 1:length(h))
  {
    tmp.term <- tmp.term + h[i] / l1.norm * ((N - (num.zeros+i) + 1/2) / N)
  }
  return(1 - 2 * tmp.term)
}

#' Spatial chaos measure.
#'
#' \code{spatial.chaos} returns the 'spatial chaos' randomness measure for
#' imaging data.
#'
#' @param im 2-D numeric matrix representing the image pixel intensities.
#' @param levels numeric (default = 30). Number of histogram bins.
#' @param morph logical (default = \code{TRUE}). Whether morphological operations should
#' be applied to the binary image.
#'
#' @return A value between 0 and 1. A value close to 1 represents a high level of
#' spatial scatteredness, a value close to 0 represents a less level of spatial
#' scatteredness. Maximum possible value is 1 - 1 / (# histogram bins)
#'
#' @references Palmer, A., Phapale, P., Chernyavsky, I., Lavigne, R., Fay,
#' D., Tarasov, A., ... & Becker, M. (2017). FDR-controlled metabolite annotation
#' for high-resolution imaging mass spectrometry. Nature methods, 14(1), 57.
#'
#' @author Paolo Inglese \email{p.inglese14@imperial.ac.uk}
#'
#' @example R/examples/sparseness.R
#'
#' @seealso \link{gini.index} \link{scatter.ratio}
#' @export
#' @importFrom imager mclosing_square as.cimg
#' @importFrom SDMTools ConnCompLabel
#'
spatial.chaos <- function(im, levels = 30, morph = TRUE)
{
  im <- im / max(im)
  num.obj <- array(0, levels)
  num.pix <- sum(im != 0)
  for (n in 1:levels)
  {
    th <- n / levels
    bw <- im > th
    if (morph)
    {
      bw <- as.matrix(mclosing_square(as.cimg(bw), 3))
    }
    lbl <- ConnCompLabel(bw)
    num.obj[n] <- length(unique(lbl[lbl != 0]))
  }
  return(1 - sum(num.obj) / (num.pix * levels))
}
