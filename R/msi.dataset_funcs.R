## .normIntensity
#' @importFrom stats median
.normIntensity <- function(x, method = "median", peak.ind = NULL, zero.offset = 0)
{
  accept.method <- c("TIC", "median", "PQN")
  if (!any(method %in% accept.method))
  {
    stop(".normIntensity: Valid methods are: ", paste0(accept.method, collapse = ", "), ".")
  }
  
  if (is.null(peak.ind)) {
    peak.ind <- c(1:ncol(x))
  }
  
  x[x == 0] <- NA
  x <- switch(method,

              "TIC" = {
                x[is.na(x)] <- 0
                if (any(x == 0)) {
                  zero.offset <- 1
                  cat('IMPORTANT!!! An offset equal to 1 is added to take into account of the zeros\n')
                  x <- x + zero.offset
                }
                cat('IMPORTANT!!! Use CLR transformation for proportional data calling varTransform(object, method = "clr")\n')
                for (i in 1:nrow(x))
                {
                  tic.value <- sum(x[i, peak.ind], na.rm = T)
                  if (tic.value == 0) {
                    stop("Error: scaling factor is 0!")
                  }
                  x[i, ] <- x[i, ] / tic.value
                }
                x
              },

              "median" = {
                for (i in 1:nrow(x))
                {
                  med.value <- median(x[i, peak.ind], na.rm = T)
                  if (is.na(med.value)) {
                    stop("Error: scaling factor is 0!")
                  }
                  x[i, ] <- x[i, ] / med.value
                }
                x
              },

              "PQN" = {
                if (all(peak.ind == c(1:ncol(x)))) {
                  warning("PQN can be used only using all peaks")
                }
                for (i in 1:nrow(x))
                {
                  tic.value <- sum(x[i, peak.ind], na.rm = T)
                  if (tic.value == 0) {
                    stop("Error: scaling factor is 0!")
                  }
                  x[i, ] <- x[i, ] / tic.value
                }
                x[x == 0] <- NA

                ## Try to use the faster method, if there is not enough RAM, then
                ## use the slower method

                ## Reference spectrum = non-zero median peaks
                ref.spectrum <- tryCatch(apply(x, 2, median, na.rm = T),
                                         error = function(e)
                                         {
                                           warning("Low memory. Using the slower
                                                   method to calculate the reference spectrum.")
                                           z <- array(NA, ncol(x))
                                           for (j in 1:ncol(x))
                                           {
                                             z[j] <- median(x[, j], na.rm = T)
                                           }
                                           return(z)
                                         })
                ## Quotients
                quotients <- tryCatch(x / rep(ref.spectrum, each = nrow(x)),
                                      error = function(e)
                                      {
                                        warning("Low memory. Using the slower
                                                method to calculate the quotients.")
                                        z <- matrix(NA, nrow(x), ncol(x))
                                        for (j in 1:nrow(x))
                                        {
                                          z[j, ] <- x[j, ] / ref.spectrum
                                        }
                                        return(z)
                                      })
                quotients[quotients == 0] <- NA
                ## Scaling factors
                sc.factor <- tryCatch(apply(quotients, 1, median, na.rm = T),
                                      error = function(e)
                                      {
                                        warning("Low memory. Using the slower
                                                method to calculate the scaling factors.")
                                        z <- array(NA, nrow(quotients))
                                        for (j in 1:nrow(quotients))
                                        {
                                          z[j] <- median(quotients[j, ], na.rm = T)
                                        }
                                        return(z)
                                      })
                rm(quotients)
                ## Normalized intensities
                x <- tryCatch({
                  sc.factor.mat <- sapply(sc.factor, function(z) rep(z, ncol(x)))
                  x / t(sc.factor.mat)
                },
                error = function(e)
                {
                  warning("Low memory. Using the slower method to calculate the
                          normalized intensities.")
                  for (j in 1:nrow(x))
                  {
                    x[j, ] <- x[j, ] / sc.factor[j]
                  }
                  return(x)
                })
                
                x[is.na(x)] <- 0
                
                x
              }
  )
  
  x[is.na(x)] <- 0
  
  return(x)
}

## Reduce heteroscedasticity
.varTransf <- function(x, method = "log")
{
  ## Check if NAs are present
  if (any(is.na(x)))
  {
    stop("NAs values found in the matrix.")
  }
  ## Check if negative values are present
  if (min(x) < 0)
  {
    stop("found negative values in the matrix.")
  }
  ## If the smallest intensity is not zero, show a warning saying that the intensities
  ## will be still summed to 1
  if (min(x) > 0 && method %in% c("log", "log2", "log10"))
  {
    warning(paste0("The smallest value is ", min(x),
                   ", however still adding 1 before transforming the values."))
  }
  
  accept.method <- c("log", "log2", "log10", "sqrt", "clr")
  if (!any(method %in% accept.method))
  {
    stop("Valid methods are:", paste0(accept.method, collapse = ", "), ".")
  }
  x <- switch(method,
              "log" = log(x + 1),
              "log2" = log2(x + 1),
              "log10" = log10(x + 1),
              "sqrt" = sqrt(x),
              "clr" = log(x) - apply(log(x), 1, mean)
             )
  return(x)
}
