## set generics for msi.dataset-class methods

if (is.null(getGeneric("getMZ")))
  setGeneric("getMZ", function(object, ...) standardGeneric("getMZ"))
if (is.null(getGeneric("getIntensityMat")))
  setGeneric("getIntensityMat", function(object, ...) standardGeneric("getIntensityMat"))
if (is.null(getGeneric("getShapeMSI")))
  setGeneric("getShapeMSI", function(object, ...) standardGeneric("getShapeMSI"))
if (is.null(getGeneric("binKmeans")))
  setGeneric("binKmeans", function(object, ...) standardGeneric("binKmeans"))
if (is.null(getGeneric("normIntensity")))
  setGeneric("normIntensity", function(object, ...) standardGeneric("normIntensity"))
if (is.null(getGeneric("varTransform")))
  setGeneric("varTransform", function(object, ...) standardGeneric("varTransform"))
if (is.null(getGeneric("applyPeaksFilter")))
  setGeneric("applyPeaksFilter", function(object, ...) standardGeneric("applyPeaksFilter"))


## Methods ---------------------------------------------------------------------

#' Return the m/z vector.
#'
#' @param object \link{msi.dataset-class} object.
#'
#' @return vector containing the m/z values.
#'
#' @example R/examples/msiDataset_getters.R
#'
#' @export
#' @aliases getMZ
#'
setMethod(f = "getMZ",
          signature = signature(object = "msi.dataset"),
          definition = function(object)
          {
            object@mz
          })

#' Return the peaks intensity matrix.
#'
#' @param object \link{msi.dataset-class} object.
#'
#' @return peaks intensity matrix. Rows represent pixels, and columns represent
#' peaks.
#'
#' @example R/examples/msiDataset_getters.R
#'
#' @export
#' @aliases getIntensityMat
#'
setMethod(f = "getIntensityMat",
          signature = signature(object = "msi.dataset"),
          definition = function(object)
          {
            object@matrix
          })

#' Return a binary mask generated applying k-means clustering
#' on peaks intensities.
#'
#' @param object \link{msi.dataset-class} object
#'
#' @return \link{ms.image-class} object representing the binary mask image.
#'
#' @example R/examples/msiDataset_binKmeans.R
#'
#' @export
#' @importFrom stats kmeans
#' @aliases binKmeans
#'
setMethod(f = "binKmeans",
          signature = signature(object = "msi.dataset"),
          definition = function(object)
          {
            y.clust <- kmeans(object@matrix, centers = 2)
            y.clust <- (y.clust$cluster == 2) * 1

            values <- matrix(y.clust, object@nrow, object@ncol)

            bw <- msImage(values, "ROI")
            bw
          }
)

#' Normalize the peaks intensities.
#'
#' @param object \link{msi.dataset-class} object.
#' @param method String (default = \code{"median"}). the normalization method to
#' be used. Valid values are: \code{"TIC"}, \code{"median"}, or \code{"PQN"}.
#' See 'Details' section.
#'
#' @details The valid values for \code{method} are:
#' \itemize{
#'   \item \code{"TIC"}: total ion current normalization assign the sum of the
#'   peaks intensities to one.
#'   \item \code{"median"}: median of spectrum intensities is scaled to one.
#'   \item \code{"PQN"}:
#'   \enumerate{
#'     \item apply \code{"TIC"} normalization
#'     \item calculate the median reference spectrum (after removing the zeros)
#'     \item calculate the quotients of peaks intensities and reference
#'     \item calculate the median of quotients for each peak (after removing the zeros)
#'     \item divide all the peak intensities by the median of quotients
#'   }
#' }
#'
#' @return object \link{msi.dataset-class} object, with normalized peaks
#' intensities.
#'
#' @author Paolo Inglese \email{p.inglese14@imperial.ac.uk}
#'
#' @references F. Dieterle, A. Ross, G. Schlotterbeck, and Hans Senn. 2006.
#' Probabilistic quotient normalization as robust method to account for dilution
#' of complex biological mixtures. Application in 1H NMR metabonomics.
#' Analytical Chemistry 78(13): 4281-4290.
#'
#' @example R/examples/msiDataset_transform.R
#'
#' @seealso \link{msi.dataset-class}
#' @export
#' @aliases normIntensity
#'
setMethod(f = "normIntensity",
          signature = signature(object = "msi.dataset"),
          definition = function(object, method = "median")
          {
            object@matrix <- .normIntensity(object@matrix, method = method)
            object
          }
)

#' Variance stabilizing transformation.
#'
#' \code{varTransform} transforms the MS intensities in order to reduce heteroscedasticity.
#'
#' @param object \link{msi.dataset-class} object. See \link{msiDataset}.
#' @param method string (default = \code{log}). Transformation method.
#' Valid values are:
#' \itemize{
#'   \item "log": log-transformation defined as log(x + 1)
#'   \item "sqrt": square-root transformation.
#' }
#'
#' @example R/examples/msiDataset_transform.R
#'
#' @return \link{msi.dataset-class} object with transformed peaks intensities.
#' @export
#' @aliases varTransform
#'
setMethod(f = "varTransform",
          signature = signature(object = "msi.dataset"),
          definition = function(object, method = "log")
          {
            object@matrix <- .varTransf(object@matrix, method)
            object
          })

#' Apply the results of a peaks filter.
#'
#' \code{applyPeaksFilter} select the peaks returned by a peak filter. Custom
#' filters can be created passing a named array of selected peak indices to
#' \link{createPeaksFilter}. Names correspond to the m/z values of the selected
#' peaks and must coincide with those of the MS dataset.
#'
#' @param object \link{msi.dataset-class} object.
#' @param peakFilter peaks filter results.
#'
#' @return \link{msi.dataset-class} object with only selected peaks.
#'
#' @example R/examples/filter_csr.R
#'
#' @aliases applyPeaksFilter-msi.dataset-method applyPeaksFilter
#' @export
#' @aliases applyPeaksFilter
#'
setMethod(f = "applyPeaksFilter",
          signature = signature(object = "msi.dataset"),
          definition = function(object, peakFilter)
          {
            .stopIfNotValidPeakFilter(peakFilter)

            # Check that the peak filter m/z values correspond to the same
            # of the dataset
            if (attr(peakFilter, "filter") == "splitPeaks")
            {
              tmp.matrix <- object@matrix
              x.merged <- matrix(NA, nrow(object@matrix), length(peakFilter$merged.peaks))

              for (i in 1:length(peakFilter$merged.peaks))
              {
                x.merged[, i] <- apply(object@matrix[, peakFilter$merged.peaks[[i]]], 1, sum)
              }
              mz.merged <- names(peakFilter$merged.peaks)
              tmp.matrix <- tmp.matrix[, -unique(unlist(peakFilter$merged.peaks, use.names = F))]
              tmp.mz <- object@mz[-unique(unlist(peakFilter$merged.peaks, use.names = F))]

              # Update with the merged values
              tmp.matrix <- cbind(tmp.matrix, x.merged)
              tmp.mz <- c(tmp.mz, mz.merged)

              # Sort the m/z values
              tmp.mz <- sort(tmp.mz, index.return = T)
              tmp.matrix <- tmp.matrix[, tmp.mz$ix]

              object@matrix <- tmp.matrix
              object@mz <- as.numeric(tmp.mz$x)
            } else
            {
              if (!any(object@mz[peakFilter$sel.peaks] %in% names(peakFilter$sel.peaks)))
              {
                stop("peakFilter not compatible with the MSI dataset.")
              }
              object@matrix <- object@matrix[, peakFilter$sel.peaks]
              object@mz <- object@mz[peakFilter$sel.peaks]
            }

            object
          }
)

#' Returns the geometrical shape of MS image
#'
#' @param object \link{msi.dataset-class} object.
#'
#' @return number of rows ans number of columns of the MS image.
#'
#' @example R/examples/msiDataset_getters.R
#'
#' @export
#' @aliases getShapeMSI
#'
setMethod(f = "getShapeMSI",
          signature = signature(object = "msi.dataset"),
          definition = function(object)
          {
            c(object@nrow, object@ncol)
          })
