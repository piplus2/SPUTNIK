#' Load the example MALDI-MSI data.
#'
#' Loads a single mouse urinary bladder MALDI mass spectrometry imaging dataset
#' acquired in positive ionization mode using Thermo qExactive Orbitrap. The
#' dataset is available at
#' "https://raw.github.com/paoloinglese/SPUTNIKexamples/master/data/bladder_maldi_prepr_MALDIquant.RData"
#' The dataset is loaded in the R environment under the variable name \code{maldiData}.
#'
#' @param verbose Logical (default = TRUE). Show additional output text.
#'
#' @return desiData MS intensity matrix. Rows represent pixels, columns represent
#' matched peaks.
#'
#' @references Rompp, A., Guenther, S., Schober, Y., Schulz, O., Takats, Z., Kummer, W.,
#' & Spengler, B. (2010). Histology by mass spectrometry: label-free tissue
#' characterization obtained from high-accuracy bioanalytical imaging. Angewandte
#' chemie international edition, 49(22), 3834-3838.
#'
#' @importFrom utils download.file
#' @export
#' @aliases bladderMALDIRompp2010
#'
bladderMALDIRompp2010 <- function(verbose = TRUE) {
  maldiData <- NULL

  if (verbose) {
    cat("Downloading the data from the repository...\n")
  }

  addr <- "https://raw.github.com/paoloinglese/SPUTNIKexamples/master/data/maldiData.rda"
  z <- paste0(tempfile(), ".rda")
  tryCatch(download.file(addr, z, mode = "wb"),
    error = function(e) {
      stop(e)
    }
  )

  if (verbose) {
    cat("Loading the data in the R environment...\n")
  }

  load(z)
  ## Delete the temporary file
  file.remove(z)
  # Return the intensity matrix
  return(maldiData)
}

#' Load the example DESI-MSI data.
#'
#' Loads a single human ovarian cancer DESI mass spectrometry imaging dataset
#' acquired in negative ionization mode using Waters XEVO-GS2 qToF. The
#' dataset is available at
#' "https://raw.github.com/paoloinglese/SPUTNIKexamples/master/data/ovarian_xevo_prepr_MALDIquant.RData"
#' The dataset is loaded in the R environment under the variable name \code{maldiData}.
#'
#' @param verbose Logical (default = TRUE). Show additional output text.
#'
#' @return maldiData MS intensity matrix. Rows represent pixels, columns represent
#' matched peaks.
#'
#' @references Doria, M. L., McKenzie, J. S., Mroz, A., Phelps, D. L., Speller,
#' A., Rosini, F., ... & Ghaem-Maghami, S. (2016). Epithelial ovarian carcinoma
#' diagnosis by desorption electrospray ionization mass spectrometry imaging.
#' Scientific reports, 6, 39219.
#'
#' @importFrom utils download.file
#' @export
#' @aliases ovarianDESIDoria2016
#'
ovarianDESIDoria2016 <- function(verbose = TRUE) {
  desiData <- NULL

  if (verbose) {
    cat("Downloading the data from the repository...\n")
  }

  addr <- "https://raw.github.com/paoloinglese/SPUTNIKexamples/master/data/desiData.rda"
  z <- paste0(tempfile(), ".RData")
  tryCatch(download.file(addr, z, mode = "wb"),
    error = function(e) {
      stop(e)
    }
  )

  if (verbose) {
    cat("Loading the data in the R environment...\n")
  }

  load(z)
  ## Delete the temporary file
  file.remove(z)
  # Return the intensity matrix
  return(desiData)
}
