## .normIntensity
#' @importFrom stats median
#' @import edgeR
.normIntensity <- function(x, method = "median", peak.ind = NULL, zero.offset = 0) {
  accept.method <- c("median", "PQN", "TIC", "TMM", "upperQuartile")
  if (!any(method %in% accept.method)) {
    stop(".normIntensity: Valid methods are: ", paste0(accept.method, collapse = ", "), ".")
  }

  if (is.null(peak.ind)) {
    peak.ind <- c(1:ncol(x))
  }

  x[x == 0] <- NA
  x <- switch(method,

    "TMM" = {
      cat('IMPORTANT!!! TMM requires log2 transformation through varTransform(object, method = "log2").\n')

      # Using the default parameters from edgeR
      logratioTrim <- 0.3
      sumTrim <- 0.05
      doWeighting <- TRUE
      Acutoff <- -1e10

      x[is.na(x)] <- 0
      x <- t(x) # TMM requires samples along columns
      nsamples <- ncol(x)

      # Check lib.size
      # Force lib.size = NULL
      lib.size <- NULL
      if (is.null(lib.size)) {
        lib.size <- colSums(x)
      } else {
        if (anyNA(lib.size)) stop("NA lib.sizes not permitted")
        if (length(lib.size) != nsamples) {
          if (length(lib.size) > 1L) {
            warning("calcNormFactors: length(lib.size) doesn't match number of samples", call. = FALSE)
          }
          lib.size <- rep(lib.size, length = nsamples)
        }
      }

      # Remove all zero rows
      allzero <- .rowSums(x > 0, nrow(x), nsamples) == 0L
      if (any(allzero)) {
        x <- x[!allzero, , drop = FALSE]
      }

      # Degenerate cases
      if (nrow(x) == 0 || nsamples == 1) {
        stop("Too many empty samples. Cannot apply TMM.")
      }

      f75 <- .calcFactorQuantile(data = x, lib.size = lib.size, p = 0.75)
      # Force refColumns = NULL
      refColumn <- NULL
      if (is.null(refColumn)) {
        refColumn <- which.min(abs(f75 - mean(f75)))
      }
      if (length(refColumn) == 0L | refColumn < 1 | refColumn > nsamples) refColumn <- 1L
      f <- rep(NA, nsamples)
      for (i in 1:nsamples) {
        f[i] <- .calcFactorTMM(
          obs = x[, i], ref = x[, refColumn], libsize.obs = lib.size[i],
          libsize.ref = lib.size[refColumn], logratioTrim = logratioTrim,
          sumTrim = sumTrim, doWeighting = doWeighting, Acutoff = Acutoff
        )
      }
      f <- f / exp(mean(log(f)))

      # Output
      names(f) <- colnames(x)

      x <- x / f
      t(x) # Restore the original order (rows = samples)
    },

    "upperQuartile" = {
      if (zero.offset != 0) {
        x[is.na(x)] <- 0
      }
      x <- x + zero.offset
      for (i in 1:nrow(x))
      {
        quartile <- quantile(x[i, ], probs = 0.75, na.rm = TRUE)
        if (is.na(quartile)) {
          warning("Upper quartile is 0!")
          quartile <- 1
        }
        x[i, ] <- x[i, ] / quartile
      }
      x
    },

    "TIC" = {
      x[is.na(x)] <- 0
      if (any(x == 0)) {
        zero.offset <- 1
        cat("IMPORTANT!!! An offset equal to 1 is added to take into account of the zeros
")
        x <- x + zero.offset
      }
      cat('IMPORTANT!!! Use CLR transformation for proportional data calling varTransform(object, method = "clr")\n')
      for (i in 1:nrow(x))
      {
        tic.value <- sum(x[i, peak.ind], na.rm = TRUE)
        if (tic.value == 0) {
          stop("TIC is 0!")
        }
        x[i, ] <- x[i, ] / tic.value
      }
      x
    },

    "median" = {
      if (zero.offset != 0) {
        x[is.na(x)] <- 0
      }
      x <- x + zero.offset
      for (i in 1:nrow(x))
      {
        med.value <- median(x[i, peak.ind], na.rm = TRUE)
        if (is.na(med.value)) {
          warning("Median is 0!")
          med.value <- 1
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
        tic.value <- sum(x[i, peak.ind], na.rm = TRUE)
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
        error = function(e) {
          warning("Low memory. Using the slower
                                                   method to calculate the reference spectrum.")
          z <- array(NA, ncol(x))
          for (j in 1:ncol(x))
          {
            z[j] <- median(x[, j], na.rm = T)
          }
          return(z)
        }
      )
      ## Quotients
      quotients <- tryCatch(x / rep(ref.spectrum, each = nrow(x)),
        error = function(e) {
          warning("Low memory. Using the slower
                                                method to calculate the quotients.")
          z <- matrix(NA, nrow(x), ncol(x))
          for (j in 1:nrow(x))
          {
            z[j, ] <- x[j, ] / ref.spectrum
          }
          return(z)
        }
      )
      quotients[quotients == 0] <- NA
      ## Scaling factors
      sc.factor <- tryCatch(apply(quotients, 1, median, na.rm = T),
        error = function(e) {
          warning("Low memory. Using the slower
                                                method to calculate the scaling factors.")
          z <- array(NA, nrow(quotients))
          for (j in 1:nrow(quotients))
          {
            z[j] <- median(quotients[j, ], na.rm = T)
          }
          return(z)
        }
      )
      rm(quotients)
      ## Normalized intensities
      x <- tryCatch({
        sc.factor.mat <- sapply(sc.factor, function(z) rep(z, ncol(x)))
        x / t(sc.factor.mat)
      },
      error = function(e) {
        warning("Low memory. Using the slower method to calculate the
                          normalized intensities.")
        for (j in 1:nrow(x))
        {
          x[j, ] <- x[j, ] / sc.factor[j]
        }
        return(x)
      }
      )

      x[is.na(x)] <- 0

      x
    }
  )

  x[is.na(x)] <- 0

  return(x)
}

## Reduce heteroscedasticity
.varTransf <- function(x, method = "log") {
  ## Check if NAs are present
  if (any(is.na(x))) {
    stop("NAs values found in the matrix.")
  }
  ## Check if negative values are present
  if (min(x) < 0) {
    stop("found negative values in the matrix.")
  }
  ## If the smallest intensity is not zero, show a warning saying that the intensities
  ## will be still summed to 1
  if (min(x) > 0 && method %in% c("log", "log2", "log10")) {
    warning(paste0(
      "The smallest value is ", min(x),
      ", however still adding 1 before transforming the values."
    ))
  }

  accept.method <- c("log", "log2", "log10", "sqrt", "clr")
  if (!any(method %in% accept.method)) {
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

##################################
## TMM normalization from edgeR ##
##################################

.calcFactorTMM <- function(obs, ref, libsize.obs = NULL, libsize.ref = NULL,
                           logratioTrim = .3, sumTrim = 0.05, doWeighting = TRUE,
                           Acutoff = -1e10)
                           # TMM between two libraries
# Mark Robinson
{
  obs <- as.numeric(obs)
  ref <- as.numeric(ref)

  if (is.null(libsize.obs)) nO <- sum(obs) else nO <- libsize.obs
  if (is.null(libsize.ref)) nR <- sum(ref) else nR <- libsize.ref

  logR <- log2((obs / nO) / (ref / nR)) # log ratio of expression, accounting for library size
  absE <- (log2(obs / nO) + log2(ref / nR)) / 2 # absolute expression
  v <- (nO - obs) / nO / obs + (nR - ref) / nR / ref # estimated asymptotic variance

  # 	remove infinite values, cutoff based on A
  fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)

  logR <- logR[fin]
  absE <- absE[fin]
  v <- v[fin]

  if (max(abs(logR)) < 1e-6) {
    return(1)
  }

  # 	taken from the original mean() function
  n <- length(logR)
  loL <- floor(n * logratioTrim) + 1
  hiL <- n + 1 - loL
  loS <- floor(n * sumTrim) + 1
  hiS <- n + 1 - loS

  # 	keep <- (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
  # 	a fix from leonardo ivan almonacid cardenas, since rank() can return
  # 	non-integer values when there are a lot of ties
  keep <- (rank(logR) >= loL & rank(logR) <= hiL) & (rank(absE) >= loS & rank(absE) <= hiS)

  if (doWeighting) {
    f <- sum(logR[keep] / v[keep], na.rm = TRUE) / sum(1 / v[keep], na.rm = TRUE)
  } else {
    f <- mean(logR[keep], na.rm = TRUE)
  }

  # 	Results will be missing if the two libraries share no features with positive counts
  # 	In this case, return unity
  if (is.na(f)) {
    f <- 0
  }
  return(2^f)
}

.calcFactorQuantile <- function(data, lib.size, p = 0.75)
                                # 	Generalized version of upper-quartile normalization
                                # 	Mark Robinson
# 	Created 16 Aug 2010
{
  # 	i <- apply(data<=0,1,all)
  # 	if(any(i)) data <- data[!i,,drop=FALSE]
  y <- t(t(data) / lib.size)
  f <- apply(y, 2, function(x) quantile(x, p = p))
  return(f)
  # 	f/exp(mean(log(f)))
}