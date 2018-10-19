## .isMSImageObject
#' @importFrom methods is
.isMSImageObject <- function(x)
{
  is(object = x, class2 = "ms.image")
}

## .isMSIDatasetObject
#' @importFrom methods is
.isMSIDatasetObject <- function(x)
{
  is(object = x, class2 = "msi.dataset")
}

## .isBinary
.isBinary <- function(x)
{
  .stopIfNotValidMSImage(x)
  if (all(unique(c(x@values)) %in% c(0, 1)))
  {
    return(TRUE)
  } else
  {
    return(FALSE)
  }
}

## .isPeakFilter
.isPeakFilter <- function(x)
{
  if (is.null(attr(x, "peak.filter")))
  {
    return(FALSE)
  }
  if (!attr(x, "peak.filter"))
  {
    return(FALSE)
  }
  return(TRUE)
}

## .stopIfNotValidMSImage
.stopIfNotValidMSImage <- function(x)
{
  if (!.isMSImageObject(x))
  {
    stop("no SPUTNIK::ms.image object!")
  }
}

## .stopIfNotValidMSIDataset
.stopIfNotValidMSIDataset <- function(x)
{
  if (!.isMSIDatasetObject(x))
  {
    stop("no SPUTNIK::msi.dataset object!")
  }
}

## .stopIfNotValidPeakFilter
.stopIfNotValidPeakFilter <- function(x)
{
  if (!.isPeakFilter(x))
  {
    stop("not a valid peak.filter result.")
  }
}

## .stopIfNotValidGlobalMethod
.stopIfNotValidGlobalMethod <- function(x)
{
  accept.values <- c("pearson", "spearman", "ssim", "nmi")
  if (!any(x %in% accept.values))
  {
    stop("Accepted method values are: ",
         paste0(accept.values, collapse = ", "), ".")
  }
}
