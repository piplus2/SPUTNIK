## set generics for msi.dataset-class methods

if (is.null(getGeneric("applyPeaksFilter"))) {
  setGeneric("applyPeaksFilter", function(object, ...) standardGeneric("applyPeaksFilter"))
}

if (is.null(getGeneric("binKmeans"))) {
  setGeneric("binKmeans", function(object, ...) standardGeneric("binKmeans"))
}

if (is.null(getGeneric("binKmeans2"))) {
  setGeneric("binKmeans2", function(object, ...) standardGeneric("binKmeans2"))
}

if (is.null(getGeneric("binSupervised"))) {
  setGeneric("binSupervised", function(object, ...) standardGeneric("binSupervised"))
}

if (is.null(getGeneric("getIntensityMat"))) {
  setGeneric("getIntensityMat", function(object, ...) standardGeneric("getIntensityMat"))
}

if (is.null(getGeneric("getMZ"))) {
  setGeneric("getMZ", function(object, ...) standardGeneric("getMZ"))
}

if (is.null(getGeneric("getShapeMSI"))) {
  setGeneric("getShapeMSI", function(object, ...) standardGeneric("getShapeMSI"))
}

if (is.null(getGeneric("normIntensity"))) {
  setGeneric("normIntensity", function(object, ...) standardGeneric("normIntensity"))
}

if (is.null(getGeneric("varTransform"))) {
  setGeneric("varTransform", function(object, ...) standardGeneric("varTransform"))
}

## Methods ---------------------------------------------------------------------

## getMZ ----

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
setMethod(
  f = "getMZ",
  signature = signature(object = "msi.dataset"),
  definition = function(object) {
    return(object@mz)
  }
)

## getIntensityMat ----

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
setMethod(
  f = "getIntensityMat",
  signature = signature(object = "msi.dataset"),
  definition = function(object) {
    x <- object@matrix
    # Remove attributes and dimnames
    for (n in names(attributes(x))[names(attributes(x)) != "dim"]) {
      attr(x, n) <- NULL
    }
    return(x)
  }
)

## binKmeans ----

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
setMethod(
  f = "binKmeans",
  signature = signature(object = "msi.dataset"),
  definition = function(object) {
    y.clust <- kmeans(object@matrix, centers = 2, iter.max = 1000, nstart = 5)
    y.clust <- (y.clust$cluster == 2) * 1

    values <- matrix(y.clust, object@nrow, object@ncol)

    return(msImage(values = values, name = "ROI", scale = FALSE))
  }
)

## binKmeans2 ----

#' Return a binary mask generated applying k-means clustering
#' on peaks intensities. A finer segmentation is obtained by using a larger
#' number of clusters than 2. The off-sample clusters are merged looking at the
#' most frequent labels in the image corners. The lookup areas are defined by
#' the kernel size.
#'
#' @param object \link{msi.dataset-class} object
#' @param mzQuery numeric. Values of m/z used to calculate the reference image.
#' 2 values are interpreted as interval, multiple or single values are searched
#' in the m/z vector. It should be left unset when using \code{useFullMZRef = TRUE}.
#' @param mzTolerance numeric (default = Inf). Tolerance in PPM to match the
#' \code{mzQueryRef}. values in the m/z vector. It overrides \code{useFullMZ}.
#' @param useFullMZ logical (default = `TRUE``). Whether all the peaks should be
#' used to calculate the reference image.
#' @param numClusters numeric (default = 4). Number of k-means clusters.
#' @param kernelSize 4D array (default = c(3, 3, 3, 3)). Array of sizes in pixels
#' of the corner kernels used to identify the off-sample clusters. The elements
#' represent the size of the top-left, top-right, bottom-right and bottom-left
#' corners. A negative value can be used to skip the corresponding corner.
#' @param numCores (default = 1). Multi-core parallel computation of k-means.
#' Each core corresponds to a repetition of k-means. If \code{numCores = 1}, a
#' serial k-means with 5 repetitions is performed.
#' @param verbose logical (default = `TRUE``). Show additional output.
#'
#' @return \link{ms.image-class} object representing the binary mask image.
#'
#' @author Paolo Inglese \email{p.inglese14@imperial.ac.uk}
#'
#' @export
#' @importFrom stats kmeans
#' @import imager parallel
#' @aliases binKmeans2
#'
setMethod(
  f = "binKmeans2",
  signature = signature(object = "msi.dataset"),
  definition = function(object,
                        mzQuery = numeric(),
                        useFullMZ = TRUE,
                        mzTolerance = Inf,
                        numClusters = 4,
                        kernelSize = c(3, 3, 3, 3),
                        numCores = 1,
                        verbose = TRUE) {
    if (length(kernelSize) == 1) {
      kernelSize <- rep(kernelSize, 4)
    }
    if (length(kernelSize) != 4) {
      stop("binKmeans2: 'kernelSize' must be a 4-elements array.")
    }
    if (any(kernelSize >= object@nrow) || any(kernelSize >= object@ncol)) {
      stop("Kernel size must be smaller than MSI shape.")
    }
    if (all(kernelSize < 1)) {
      stop("binKmeans2: at least one positive value is required for 'kernelSize'.")
    }
    if (length(mzQuery) != 0 && length(mzTolerance) == 0) {
      stop("binKmeans2: 'mzTolerance' missing.")
    }
    if (length(mzQuery) != 0 && useFullMZ) {
      stop("mzQuery and useFullMZ are incompatible.")
    }
    
    # Match the peaks indices
    if (length(mzQuery) != 0) {
      cat("Matching the query M/Z values...\n")
      mz.indices <- .mzQueryIndices(mzQuery, object@mz, mzTolerance,
                                    verbose = verbose)
    } else if (useFullMZ) {
      mz.indices <- seq(1, length(object@mz))
    } else {
      stop("Must provide mzQuery or set useFullMZ = TRUE")
    }

    # Parallel
    if (numCores > 1) {
      closeAllConnections()
      gc(reset = TRUE)
      cl <- makeCluster(numCores)

      clusterExport(
        cl = cl, varlist = c("object", "numClusters", "mz.indices"),
        envir = environment()
      )
      results <- clusterApply(
        cl, 1:10,
        function(n, x) {
          set.seed(NULL)
          kmeans(x, numClusters,
            nstart = 1,
            iter.max = 1000
          )
        },
        object@matrix[, mz.indices]
      )
      stopCluster(cl = cl)
      closeAllConnections()
      gc()
      # Get the clusters with the smallest WSS
      i <- sapply(results, function(result) result$tot.withinss)
      y.clust <- results[[which.min(i)]]
    } else
    # Serial
    {
      y.clust <- kmeans(object@matrix[, mz.indices],
        centers = numClusters,
        iter.max = 1000, nstart = 5
      )
    }

    # Merge the sample-related clusters
    
    roi <- matrix(0, object@nrow, object@ncol)
    for (k in 1:numClusters)
    {
      curr_clust_im <- matrix(0, object@nrow, object@ncol)
      curr_clust_im[y.clust$cluster == k] <- 1
      # Top-left corner
      if (kernelSize[1] > 0) {
        if (.mode(curr_clust_im[1:kernelSize[1], 1:kernelSize[1]] == 1)) {
          next()
        }
      }
      # Top-right
      if (kernelSize[2] > 0) {
        if (.mode(curr_clust_im[
          1:kernelSize[2],
          (object@ncol - kernelSize[2]):object@ncol
          ] == 1)) {
          next()
        }
      }
      # Bottom-right
      if (kernelSize[3] > 0) {
        if (.mode(curr_clust_im[
          (object@nrow - kernelSize[3]):object@nrow,
          (object@ncol - kernelSize[3]):object@ncol
        ] == 1)) {
          next()
        }
      }
      # Bottom-left
      if (kernelSize[4] > 0) {
        if (.mode(curr_clust_im[
          (object@nrow - kernelSize[4]):object@nrow,
          1:kernelSize[4]
        ] == 1)) {
          next()
        }
      }
      roi[curr_clust_im == 1] <- roi[curr_clust_im == 1] + 1
    }

    return(msImage(values = roi, name = "ROI", scale = FALSE))
  }
)

## binSupervised ----

#' Return a binary mask generated applying a supervised classifier
#' on peaks intensities using manually selected regions corresponding to off-sample
#' and sample-related areas.
#'
#' @param object \link{msi.dataset-class} object
#' @param refImage \link{ms.image-class} object. Image used as reference to
#' manually select the ROI pixels.
#' @param mzQuery numeric. Values of m/z used to calculate the reference image.
#' 2 values are interpreted as interval, multiple or single values are searched
#' in the m/z vector. It overrides \code{useFullMZ}.
#' @param mzTolerance numeric (default = Inf). Tolerance in PPM to match the
#' \code{mzQueryRef}. values in the m/z vector. It overrides \code{useFullMZ}.
#' @param useFullMZ logical (default = `TRUE``). Whether all the peaks should be
#' used to perform the supervised segmentation.
#' @param method string (default = 'svm'). Supervised classifier used to segment
#' the ROI.
#' @param verbose boolean (default = `TRUE`). Additional output.
#'
#' @return \link{ms.image-class} object representing the binary mask image.
#'
#' @author Paolo Inglese \email{p.inglese14@imperial.ac.uk}
#'
#' @export
#' @import imager e1071
#' @importFrom stats predict quantile
#' @aliases binSupervised
#'
setMethod(
  f = "binSupervised",
  signature = signature(object = "msi.dataset"),
  definition = function(object,
                        refImage,
                        mzQuery = numeric(), # Filter m/z values
                        mzTolerance = Inf, #
                        useFullMZ = TRUE, #
                        method = "svm",
                        verbose = TRUE) {
    accept.methods <- c("svm")

    .stopIfNotValidMSIDataset(object)
    .stopIfNotValidMSImage(refImage)
    stopifnot(is.numeric(mzQuery))
    stopifnot(is.logical(useFullMZ))
    stopifnot(is.numeric(mzTolerance))
    stopifnot(method %in% accept.methods)

    if (length(mzQuery) != 0 && useFullMZ) {
      stop("mzQuery and useFullMZ are incompatible.")
    }
    
    # Match the peaks indices
    if (length(mzQuery) != 0) {
      cat("Matching the query M/Z values...\n")
      mz.indices <- .mzQueryIndices(mzQuery, object@mz, mzTolerance, verbose = verbose)
    } else if (useFullMZ) {
      mz.indices <- seq(1, length(object@mz))
    } else {
      stop("Must provide mzQuery or set useFullMZ = TRUE")
    }

    # User-defined pixels
    userCoords <- vector(mode = "list", length = 2)
    names(userCoords) <- c("off-sample", "sample")

    for (i in 1:length(userCoords))
    {
      cat(paste0("Select the ", names(userCoords)[i], " area..."))
      userCoords[[i]] <- grabRect(as.cimg(refImage@values), output = "coord")
    }

    # Define the mask corresponding to the user-defined pixels
    mask <- matrix(0, object@nrow, object@ncol)
    mask[
      seq(userCoords[[1]][1], userCoords[[1]][3]),
      seq(userCoords[[1]][2], userCoords[[1]][4])
    ] <- 1
    mask[
      seq(userCoords[[2]][1], userCoords[[2]][3]),
      seq(userCoords[[2]][2], userCoords[[2]][4])
    ] <- 2

    # Classify the pixels
    idx.train <- which(mask != 0)
    idx.test <- which(mask == 0)

    y <- factor(mask[idx.train])
    stopifnot(all(sort(unique(y)) == c(1, 2)))
    
    cat("Segmentation...")
    mdl <- switch(method,
      "svm" = svm(object@matrix[idx.train, ], y, kernel = "linear")
    )

    ypred <- as.character(c(mask))
    ypred[idx.test] <- as.character(predict(mdl, object@matrix[idx.test, ]))

    ypred <- as.numeric(ypred)
    stopifnot(all(sort(unique(ypred)) == c(1, 2)))

    binRoi <- (ypred == 2) * 1
    binRoi <- matrix(binRoi, object@nrow, object@ncol)

    return(msImage(values = binRoi, name = "ROI", scale = FALSE))
  }
)

## normIntensity ----

#' Normalize the peaks intensities.
#'
#' @param object \link{msi.dataset-class} object.
#' @param method string (default = \code{"median"}). The normalization method to
#' be used. Valid values are: \code{"median"}, \code{"PQN"}, \code{"TIC"},
#' \code{TMM}, or \code{"upperQuartile"}.
#' See 'Details' section.
#' @param peaksInd numeric array (default = NULL). Array of peak indices used to
#' calculate the scaling factors (TIC, median). If NULL, all the peaks are used.
#' @param offsetZero numeric (default = 0). This value is added to all the peak
#' intensities to take into accounts of the zeros.
#'
#' @details The valid values for \code{method} are:
#' \itemize{
#'   \item \code{"median"}: median of spectrum intensities is scaled to one.
#'   \item \code{"PQN"}:
#'   \enumerate{
#'     \item apply \code{"TIC"} normalization
#'     \item calculate the median reference spectrum (after removing the zeros)
#'     \item calculate the quotients of peaks intensities and reference
#'     \item calculate the median of quotients for each peak (after removing the zeros)
#'     \item divide all the peak intensities by the median of quotients
#'   }
#'   \item \code{"TIC"}: total ion current normalization assign the sum of the
#'   peaks intensities to one.
#'   \item \code{"TMM"}: trimmed mean of M-values (TMM with zero pairing).
#'   Called TMMwzp in edgeR.
#'   \item \code{"upperQuartile"}: spectra are scaled by their 3rd quartile.
#' }
#'
#' @return object \link{msi.dataset-class} object, with normalized peaks
#' intensities.
#'
#' When using TIC scaling, if zeros are present in the matrix, a positive offset
#' must be added to all the peak intensities through the parameter \code{offsetZero}.
#' This is necessary for applying the CLR transformation. TIC scaling transforms the
#' spectra into compositional data; in this case the CLR transformation must be
#' applied through the varTransform function.
#'
#' @author Paolo Inglese \email{p.inglese14@imperial.ac.uk}
#'
#' @references F. Dieterle, A. Ross, G. Schlotterbeck, and Hans Senn. 2006.
#' Probabilistic quotient normalization as robust method to account for dilution
#' of complex biological mixtures. Application in 1H NMR metabonomics.
#' Analytical Chemistry 78(13): 4281-4290.
#' @references Robinson MD, Oshlack A (2010). A scaling normalization method for
#' differential expression analysis of RNA-seq data. Genome Biology 11, R25.
#'
#' @example R/examples/msiDataset_transform.R
#'
#' @seealso \link{msi.dataset-class}
#' @export
#' @aliases normIntensity
#'
setMethod(
  f = "normIntensity",
  signature = signature(object = "msi.dataset"),
  definition = function(object, method = "median", peaksInd = NULL, offsetZero = 0) {
    
    if (object@norm != "none") {
      stop("MSI already normalized. Create a new msiDataset object to use a
           different normalization method.")
    }
    if (length(offsetZero) > 1 || !is.numeric(offsetZero)) {
      stop("offsetZero must be a numeric value.")
    }
    if (is.na(offsetZero) || is.infinite(offsetZero)) {
      stop("offsetZero must be finite.")
    }
    
    object@matrix <- .normIntensity(object@matrix,
      method = method,
      peak.ind = peaksInd,
      zero.offset = offsetZero
    )
    object@norm <- method
    object@normoffset <- offsetZero

    return(object)
  }
)

## varTransform ----

#' Variance stabilizing transformation.
#'
#' \code{varTransform} transforms the MS intensities in order to reduce heteroscedasticity.
#'
#' @param object \link{msi.dataset-class} object. See \link{msiDataset}.
#' @param offsetZero numeric (default = 1). This value is added to all the peak
#' intensities to take into accounts of the zeros. It must be positive.
#' @param method string (default = \code{log}). Transformation method.
#' Valid values are:
#' \itemize{
#'   \item "log", "log2", "log10": log-transformation defined as \code{log(x + offsetZero)}.
#'   \item "sqrt": square-root transformation.
#'   \item "clr": centered log-transformation. To be used when TIC scaling
#'   normalization is applied.
#' }
#'
#' @example R/examples/msiDataset_transform.R
#'
#' @return \link{msi.dataset-class} object with transformed peaks intensities.
#' @export
#' @aliases varTransform
#'
setMethod(
  f = "varTransform",
  signature = signature(object = "msi.dataset"),
  definition = function(object, method = "log", offsetZero = 1) {
    
    if (length(offsetZero) > 1 || !is.numeric(offsetZero)) {
      stop("offsetZero must be a number.")
    }
    if (is.na(offsetZero) || is.infinite(offsetZero)) {
      stop("offsetZero must be finite.")
    }
    if (offsetZero < 0) {
      stop("offsetZero must be positive.")
    }
    
    if (object@vartr != "none") {
      stop("MSI already transformed. Create a new msiDataset object to use a 
           different transformation method.")
    }
    object@matrix <- .varTransf(object@matrix,
                                method = method,
                                zero.offset = offsetZero,
                                norm.method = object@norm)
    object@vartr <- method
    object@vartroffset <- offsetZero

    return(object)
  }
)

## applyPeaksFilter ----

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
setMethod(
  f = "applyPeaksFilter",
  signature = signature(object = "msi.dataset"),
  definition = function(object, peakFilter) {
    .stopIfNotValidPeakFilter(peakFilter)

    # Check that the peak filter m/z values correspond to the same
    # of the dataset
    if (attr(peakFilter, "filter") == "splitPeaks") {
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
    } else {
      if (!any(object@mz[peakFilter$sel.peaks] %in% names(peakFilter$sel.peaks))) {
        stop("peakFilter is not compatible with the MSI dataset.")
      }
      object@matrix <- object@matrix[, peakFilter$sel.peaks]
      object@mz <- object@mz[peakFilter$sel.peaks]
    }

    return(object)
  }
)

## getShapeMSI ----

#' Returns the geometrical shape of MSI dataset
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
setMethod(
  f = "getShapeMSI",
  signature = signature(object = "msi.dataset"),
  definition = function(object) {
    image.shape <- c(object@nrow, object@ncol)
    names(image.shape) <- c("nrow", "ncol")
    return(image.shape)
  }
)
