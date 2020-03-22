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
  peak.ind <- sort(unique(peak.ind))

  if (zero.offset != 0) {
    x[is.na(x)] <- 0
    x <- x + zero.offset
  }
  
  x[x == 0] <- NA
  x <- switch(method,

    "TMM" = {
      # Using the default parameters from edgeR
      logratioTrim <- 0.3
      sumTrim <- 0.05
      doWeighting <- TRUE
      Acutoff <- -1e10

      x[is.na(x)] <- 0
      
      # Remove empty pixels
      empty.pixel <- rowSums(x) == 0
      if (sum(empty.pixel)) {
        warning(cat(sum(empty.pixel), " empty pixels found.\n"))
        # Save the original number of pixels for final reconstruction
        orig.num.pixels <- nrow(x)
        x <- x[!empty.pixel, ]
      }
      
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
        f[i] <- .calcFactorTMMwsp(
          obs = x[, i], ref = x[, refColumn], libsize.obs = lib.size[i],
          libsize.ref = lib.size[refColumn], logratioTrim = logratioTrim,
          sumTrim = sumTrim, doWeighting = doWeighting, Acutoff = Acutoff
        )
      }
      f <- f / exp(mean(log(f)))

      # Output
      names(f) <- colnames(x)

      x <- x / f
      x <- t(x) # Restore the original order (rows = samples)
      
      # If empty pixels, restore the original dimensions
      if (sum(empty.pixel) > 0) {
        X <- matrix(0, orig.num.pixels, ncol(x))
        X[!empty.pixel, ] <- x
        x <- X
        rm(X)
        gc()
      }
      
      x
    },

    "upperQuartile" = {
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
      if (any(is.na(x)) && zero.offset == 0) {
        warning("Found peaks with zero intensity. Remember to add a numeric offset to the variance stabilizing transformation.")
      }
      cat('IMPORTANT!!! Use CLR transformation for proportional data calling varTransform(object, method = "clr")\n')
      for (i in 1:nrow(x))
      {
        tic.value <- sum(x[i, peak.ind], na.rm = TRUE)
        if (tic.value == 0) {
          warning(paste0("Pixel ", i, "has TIC = 0. The spectrum will not be scaled"))
        }
        x[i, ] <- x[i, ] / tic.value
      }
      x
    },

    "median" = {
      for (i in 1:nrow(x))
      {
        med.value <- median(x[i, peak.ind], na.rm = TRUE)
        if (is.na(med.value)) {
          warning(paste0("Pixel ", i, " has median intensiy equal to 0. The spectrum will not be scaled."))
          med.value <- 1
        }
        x[i, ] <- x[i, ] / med.value
      }
      x
    },

    "PQN" = {
      if (!all(peak.ind == c(1:ncol(x)))) {
        warning("PQN can be only applied on all peaks")
      }
      for (i in 1:nrow(x))
      {
        tic.value <- sum(x[i, ], na.rm = TRUE)
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
          warning("Low memory. Using the slower method to calculate the reference spectrum.")
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
          warning("Low memory. Using the slower method to calculate the quotients.")
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
          warning("Low memory. Using the slower method to calculate the scaling factors.")
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
        warning("Low memory. Using the slower method to calculate the normalized intensities.")
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

  # Setting all missing values to 0
  x[is.na(x)] <- 0

  return(x)
}

## Reduce heteroscedasticity
.varTransf <- function(x, method = "log", zero.offset = 0, norm.method) {
  ## Check if NAs are present
  if (any(is.na(x))) {
    stop("NAs values found in the matrix.")
  }
  ## Check if negative values are present
  if (min(x) < 0) {
    stop("found negative values in the matrix.")
  }
  # Check the selected method
  accept.method <- c("log", "log2", "log10", "sqrt", "clr")
  if (!any(method %in% accept.method)) {
    stop("Valid methods are:", paste0(accept.method, collapse = ", "), ".")
  }
  
  x <- x + zero.offset
  
  ## If the smallest intensity is not zero, show a warning saying that the intensities
  ## will be still summed to 1
  if (min(x) == 0 && method %in% c("log", "log2", "log10", "clr")) {
    stop("Method ", method, " cannot be applied if zeros are present. Add an numeric offset.")
  }
  
  # CLR must be used for TIC normalization
  if (norm.method == "TIC" && method != "clr") {
    stop("'clr' must be used with TIC scaling normalized data!!!")
  }
  
  x <- switch(method,
    "log" = log(x),
    "log2" = log2(x),
    "log10" = log10(x),
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

.calcFactorTMMwsp <- function(obs, ref, libsize.obs = NULL, libsize.ref = NULL, logratioTrim = .3, sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e10)
                              # 	TMM with pairing of singleton positive counts between the obs and ref libraries
                              # 	Gordon Smyth
# 	Created 19 Sep 2018. Last modified 23 April 2019.
{
  obs <- as.numeric(obs)
  ref <- as.numeric(ref)

  # 	epsilon serves as floating-point zero
  eps <- 1e-14

  # 	Identify zero counts
  pos.obs <- (obs > eps)
  pos.ref <- (ref > eps)
  npos <- 2L * pos.obs + pos.ref

  # 	Remove double zeros and NAs
  i <- which(npos == 0L | is.na(npos))
  if (length(i)) {
    obs <- obs[-i]
    ref <- ref[-i]
    npos <- npos[-i]
  }

  # 	Check library sizes
  if (is.null(libsize.obs)) libsize.obs <- sum(obs)
  if (is.null(libsize.ref)) libsize.ref <- sum(ref)

  # 	Pair up as many singleton positives as possible
  # 	The unpaired singleton positives are discarded so that no zeros remain
  zero.obs <- (npos == 1L)
  zero.ref <- (npos == 2L)
  k <- (zero.obs | zero.ref)
  n.eligible.singles <- min(sum(zero.obs), sum(zero.ref))
  if (n.eligible.singles > 0L) {
    refk <- sort(ref[k], decreasing = TRUE)[1:n.eligible.singles]
    obsk <- sort(obs[k], decreasing = TRUE)[1:n.eligible.singles]
    obs <- c(obs[!k], obsk)
    ref <- c(ref[!k], refk)
  } else {
    obs <- obs[!k]
    ref <- ref[!k]
  }

  # 	Any left?
  n <- length(obs)
  if (n == 0L) {
    return(1)
  }

  # 	Compute M and A values
  obs.p <- obs / libsize.obs
  ref.p <- ref / libsize.ref
  M <- log2(obs.p / ref.p)
  A <- 0.5 * log2(obs.p * ref.p)

  # 	If M all zero, return 1
  if (max(abs(M)) < 1e-6) {
    return(1)
  }

  # 	M order, breaking ties by shrunk M
  obs.p.shrunk <- (obs + 0.5) / (libsize.obs + 0.5)
  ref.p.shrunk <- (ref + 0.5) / (libsize.ref + 0.5)
  M.shrunk <- log2(obs.p.shrunk / ref.p.shrunk)
  o.M <- order(M, M.shrunk)

  # 	A order
  o.A <- order(A)

  # 	Trim
  loM <- as.integer(n * logratioTrim) + 1L
  hiM <- n + 1L - loM
  keep.M <- rep.int(FALSE, n)
  keep.M[o.M[loM:hiM]] <- TRUE
  loA <- as.integer(n * sumTrim) + 1L
  hiA <- n + 1L - loA
  keep.A <- rep.int(FALSE, n)
  keep.A[o.A[loA:hiA]] <- TRUE
  keep <- keep.M & keep.A
  M <- M[keep]

  # 	Average the M values
  if (doWeighting) {
    obs.p <- obs.p[keep]
    ref.p <- ref.p[keep]
    v <- (1 - obs.p) / obs.p / libsize.obs + (1 - ref.p) / ref.p / libsize.ref
    w <- (1 + 1e-6) / (v + 1e-6)
    TMM <- sum(w * M) / sum(w)
  } else {
    TMM <- mean(M)
  }

  2^TMM
}
