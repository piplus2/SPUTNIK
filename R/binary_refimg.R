#' Calculate the binary reference image using Otsu's thresholding.
#'
#' @param image \link{ms.image-class} object. See \link{msImage}.
#'
#' @return \link{ms.image-class} object with binary intensities.
#'#'
#' @export
#' @aliases refImageBinaryOtsu
#'
refImageBinaryOtsu <- function(image) {
  .stopIfNotValidMSImage(image)
  return(binOtsu(image))
}


#' Calculate the binary reference image using k-means clustering. K-Means is run
#' on the first `npcs` principal components to speed up the calculations.
#'
#' @param dataset \link{msi.dataset-class} object. See \link{msiDataset}.
#' @param npcs int (default = 10). Number of principal components to calculate.
#' @param alignTo string (default = "detected"). Sample reference image to align
#' the estimate binary image. It is expected to correlate with the sample location.
#' @param invertAligned boolean (default = FALSE). If TRUE, the binary image is
#' inverted after being aligned to the sample reference (\code{alignTo}).
#'
#' @return \link{ms.image-class} object with binary intensities.
#'
#' @export
#' @aliases refImageOtsu
#'
refImageBinaryKmeans <- function(dataset, npcs = 10, alignTo = "detected",
                                 invertAligned = FALSE) {
  .stopIfNotValidMSIDataset(dataset)
  return(binKmeans(dataset, npcs = npcs, ref = alignTo, invert = invertAligned))
}


#' Calculate the binary reference image using k-means clustering with multi-cluster
#' merging. K-means is run on the first `npcs` principal 
#' components to speed up the calculations.
#'
#' @param dataset \link{msi.dataset-class} object. See \link{msiDataset}.
#' @param npcs int (default = 10). Number of principal components to calculate.
#' @param mzQuery numeric. Values of m/z used to calculate the reference image.
#' 2 values are interpreted as interval, multiple or single values are searched
#' in the m/z vector. It overrides the argument \code{useFullMZ}.
#' @param mzTolerance numeric (default = Inf). Tolerance in PPM to match the
#' \code{mzQueryRef} values in the m/z vector.
#' @param useFullMZ logical (default = TRUE). Whether all the peaks should be
#' used to calculate the reference image.
#' @param numClusters numeric (default = 4). Number of clusters.
#' @param kernelSize 4-D numeric array or numeric (default = 5).
#' Each element of the 4-D array represents the size of the corners square kernels
#' used to determine the off-tissue clusters. The element order is clockwise:
#' top-left, top-right, bottom-left, bottom-right. If negative, the corresponding
#' corner is skipped. If only a single value is passed, the same kernel size is
#' used for the 4 corners.
#' @param cores numeric (default = 1). Number of CPU cores for parallel k-means.
#' @param verbose boolean (default = TRUE). Additional output.
#'
#' @return \link{ms.image-class} object with binary intensities.
#'
#' @export
#' @aliases refImageBinaryKmeansMulti
#'
refImageBinaryKmeansMulti <- function(dataset,
                                      npcs = 10,
                                      mzQuery = numeric(),
                                      mzTolerance = Inf,
                                      useFullMZ = TRUE, 
                                      numClusters = 4,
                                      kernelSize = 5,
                                      cores = 1,
                                      verbose = TRUE) {
  .stopIfNotValidMSIDataset(dataset)
  
  bw <- binKmeans2(dataset, mzQuery = mzQuery, useFullMZ = useFullMZ,
                   mzTolerance = mzTolerance, numClusters = numClusters,
                   kernelSize = kernelSize, numCores = cores, verbose = verbose)
  
  return (bw)
}


#' Calculate the binary reference image using linear SVM trained on manually 
#' selected pixels.
#'
#' @param dataset \link{msi.dataset-class} object. See \link{msiDataset}.
#' @param mzQueryRef numeric. Values of m/z used to calculate the reference image.
#' 2 values are interpreted as interval, multiple or single values are searched
#' in the m/z vector. It overrides the argument \code{useFullMZ}.
#' @param mzTolerance numeric (default = Inf). Tolerance in PPM to match the
#' \code{mzQueryRef} values in the m/z vector.
#' @param useFullMZ logical (default = TRUE). Whether all the peaks should be
#' used to calculate the reference image.
#' 
#' @return \link{ms.image-class} object with binary intensities.
#'
#' @import irlba
#' @export
#' @aliases refImageBinarySVM

refImageBinarySVM <- function(dataset,
                              mzQueryRef = numeric(),
                              mzTolerance = Inf,
                              useFullMZ = TRUE) {
  pc <- irlba::prcomp_irlba(dataset@matrix, center = TRUE, scale. = TRUE, n = 1)
  ref.img <- msImage(matrix(pc$x[, 1], dataset@nrow, dataset@ncol),
                     name = 'PC1', scale = TRUE)
  bw <- binSupervised(dataset,
                      refImage = ref.img,
                      mzQuery = mzQueryRef,
                      useFullMZ = useFullMZ,
                      mzTolerance = mzTolerance,
                      method = "svm")
}
