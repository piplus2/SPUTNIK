
#' \code{refImageContinuous} returns the reference image, calculated using the
#' \code{method}.
#' This image represents the basic measure for the filters in SPUTNIK.
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
#' @param mzQueryRef numeric. Mass-to-charge ratios used to calculate the reference image.
#' Two values are interpreted as interval, multiple or single values are searched
#' in the m/z vector. It overrides the \code{useFullMZRef} argument.
#' @param mzTolerance numeric (default = Inf). Search window in parts-per-million units
#' for the \code{mzQueryRef}.
#' @param useFullMZRef logical (default = \code{TRUE}). Whether all peaks should be
#' used to calculate the reference image. Ignored if \code{mzQueryRef} is provided.
#' @param doSmooth logical (default = FALSE). If \code{TRUE}, the reference image is
#' smoothed using a Gaussian kernel.
#' @param smoothSigma numeric (default = 2). Standard deviation of the smoothing 
#' Gaussian kernel.
#' @param alignTo string (default = "detected"). The reference image is aligned
#'to the image representing:
#' \itemize{
#'  \item "detected": number of detected peaks
#'  \item "tic": total-ion-count image
#' }
#' The reference image will have a positive Pearson's correlation with the selected image.
#' @param invertAligned logical (default = FALSE). If \code{TRUE}, the reference image
#' has negative correlation with the selected image in \code{alignTo}.
#' @param verbose logical (default = TRUE). Additional output text.
#'
#' @details Function to extract the continuous reference image from a 
#' \code{\link{msi.dataset-class}} object.
#' The continuous reference image represents the spatial location of the sample.
#' By default, it is aligned with either the image representing the number of detected
#' peaks, or the total-ion-count in all pixels. It is expected to be higher in the
#' region occupied by the sample (positive correlation with the mask representing
#' the real sample pixels). In some cases, the alignment images can have higher
#' values in the pixels outside of the sample. In these circumstances, the argument
#' \code{invertAligned} should be set to \code{TRUE}.
#'
#' @example R/examples/graph_funcs.R
#'
#' @seealso msiDataset
#' @export
#' 
refImageContinuous <- function(msiData,
                               method = "sum",
                               mzQueryRef = numeric(),
                               mzTolerance = Inf,
                               useFullMZRef = TRUE,
                               doSmooth = FALSE,
                               smoothSigma = 2,
                               alignTo = "detected",
                               invertAligned = FALSE,
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
    sampleReference = alignTo,
    invertAlign = invertAligned
  )
  
  return(ref.image)
}
