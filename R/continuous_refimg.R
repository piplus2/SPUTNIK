
#' \code{RefImageContinuous} returns the reference image, calculated using the
#' \code{refMethod}.
#' This images represents the basic measure for the filters in SPUTNIK.
#'
#' @param msiData \link[SPUTNIK]{msiDataset} object..
#' @param method string (default = "sum"). Method used to calculate the
#' reference image. Valid values are:
#' \itemize{
#'   \item "sum": peak intensities sum
#'   \item "mean": average peak intensities (without zeros)
#'   \item "median": median peak intensities (without zeros)
#'   \item "pca": first principal component scores.
#' }
#' @param mzQueryRef numeric. Values of m/z used to calculate the reference image.
#' 2 values are interpreted as interval, multiple or single values are searched
#' in the m/z vector. It overrides the param \code{useFullMZRef}.
#' @param mzTolerance numeric (default = Inf). Tolerance in PPM to match the
#' \code{mzQueryRef} values in the m/z vector.
#' @param useFullMZRef logical (default = TRUE). Whether all the peaks should be
#' used to calculate the reference image.
#' @param doSmooth logical (default = FALSE). Whether the reference image
#' should be smoothed before binarizing. Only valid for \code{roiMethod = "otsu"}.
#' @param smoothSigma numeric (default = 2). Standard deviation of Gaussian
#' kernel.
#' @param sampleReference string (default = "detected"). Image used as reference for
#' aligning the estimated reference and ROI image. Valid values are:
#' \itemize{
#'  \item "detected": number of detected peaks
#'  \item "tic": total-ion-count image
#' }
#' @param alignToReference logical (default = TRUE). Whether the reference image
#' colors should correlate with the number of detected ions. In case the number
#' of detected ions is larger outside of the sample, set it to FALSE.
#' @param verbose logical (default = TRUE). Additional output text.
#'
#' @details Function to extract the continuous reference image from a 
#' \code{\link{msi.dataset-class}} object.
#' Multiple methods can be used to extract the reference image, which afterwards 
#' can be used as argument for the \code{\link{globalPeaksFilter}} filter.
#'
#' @author Paolo Inglese \email{p.inglese14@imperial.ac.uk}
#'
#' @example R/examples/graph_funcs.R
#'
#' @seealso msiDataset
#' @export
#' 
RefImageContinuous <- function(msiData,
                            method = "sum",
                            mzQueryRef = numeric(),
                            mzTolerance = Inf,
                            useFullMZRef = TRUE,
                            doSmooth = FALSE,
                            smoothSigma = 2,
                            sampleReference = "detected",
                            alignToReference = TRUE,
                            verbose = TRUE) {
  .stopIfNotValidMSIDataset(msiData)
  
  # Reference image
  ref.image <- .refImage(
    msiData = msiData,
    method = method,
    mzQuery = mzQueryRef,
    mzTolerance = mzTolerance,
    useFullMZ = useFullMZRef,
    smoothIm = doSmooth,
    smoothSigma = smoothSigma,
    sampleReference = sampleReference,
    alignToReference = alignToReference
  )
  
  return(ref.image)
}
