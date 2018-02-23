## .normIntensity
#' @importFrom stats median
.normIntensity <- function(x, method = "median")
{
  accept.method <- c("TIC", "median", "PQN")
  if (!any(method %in% accept.method))
  {
    stop("Valid methods are: ", paste0(accept.method, collapse = ", "), ".")
  }
  x[x == 0] <- NA
  x <- switch(method,

              "TIC" = {
                for (i in 1:nrow(x))
                {
                  x[i, ] <- x[i, ] / sum(x[i, ], na.rm = T)
                }
                x
              },

              "median" = {
                for (i in 1:nrow(x))
                {
                  x[i, ] <- x[i, ] / median(x[i, ], na.rm = T)
                }
                x
              },

              "PQN" = {

                for (i in 1:nrow(x))
                {
                  x[i, ] <- x[i, ] / sum(x[i, ], na.rm = T)
                }
                x[x == 0] <- NA

                ## Try to use the faster method, if memory is not enough then
                ## use the slower method

                ## Reference spectrum = non-zero median peaks
                ref.spectrum <- tryCatch(apply(x, 2, median, na.rm = T),
                                         error = function(e) {
                                           warning("Low memory. Using the slower method to calculate the reference spectrum.")
                                           z <- array(NA, ncol(x))
                                           for (j in 1:ncol(x))
                                           {
                                             z[j] <- median(x[, j], na.rm = T)
                                           }
                                           return(z)
                                         })
                ## Quotients
                quotients <- tryCatch(x / rep(ref.spectrum, each = nrow(x)),
                                      error = function(e) {
                                        warning("Low memory. Using the slower method to calculate the quotients.")
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
                                      error = function(e) {
                                        warning("Low memory. Using the slower method to calculate the scaling factors.")
                                        z <- array(NA, nrow(quotients))
                                        for (j in 1:nrow(quotients))
                                        {
                                          z[j] <- median(quotients[j, ], na.rm = T)
                                        }
                                        return(z)
                                      })
                rm(quotients)
                ## Normalized intensities
                x <- tryCatch(x / rep(sc.factor, each = ncol(x)),
                              error = function(e) {
                                warning("Low memory. Using the slower method to calculate the normalized intensities.")
                                for (j in 1:nrow(x))
                                {
                                  x[j, ] <- x[j, ] / sc.factor[j]
                                }
                                return(x)
                              })
                x[is.na(x)] <- 0
                x
              })
  x[is.na(x)] <- 0
  x
}

## .varTransf
.varTransf <- function(x, method = "log")
{
  accept.method <- c("log", "sqrt")
  if (!any(method %in% accept.method))
  {
    stop("Valid methods are:", paste0(accept.method, collapse = ", "), ".")
  }
  x <- switch(method,
              "log" = log(x + 1),
              "sqrt" = sqrt(x))
}
