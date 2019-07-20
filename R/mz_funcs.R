## .distance.ppm
.distance.ppm <- function(mz1, mz2) {
  return(abs(mz1 - mz2) / mean(mz1, mz2) * 1e6)
}

## .distance.ppm.interval
.distance.ppm.interval <- function(mz.vec) {
  return(.distance.ppm(min(mz.vec), max(mz.vec)))
}

## .mzQueryIndices
.mzQueryIndices <- function(mz.query, mz.vector, mz.tol, verbose = FALSE) {
  mz.query <- unique(mz.query)
  search.type <- .mz.search.type(mz.query, verbose)

  mz.matched.indices <- switch(search.type,

    "interval" = {
      mz.matched <- .match.mz.interval(
        mz.query, mz.vector,
        mz.tol
      )
      if (length(unique(mz.matched)) == 1) {
        warning("Matched m/z values are identical. Using only 1 value.")
        mz.matched[1]
      } else {
        seq(mz.matched[1], mz.matched[2])
      }
    },

    "single" = {
      .match.mz.single(mz.query, mz.vector, mz.tol)
    },

    "array" = {
      mz.matched <- .match.mz.array(mz.query, mz.vector, mz.tol)
      mz.matched <- unique(mz.matched)
      if (length(mz.matched) == 1) {
        warning("Matched m/z values are identical. Using only 1 value.")
        mz.matched
      } else {
        mz.matched
      }
    }
  )

  return(mz.matched.indices)
}

## .mz.search.type
.mz.search.type <- function(mz = NULL,
                            verbose = FALSE) {
  out <- NULL

  if (is.numeric(mz) && length(unique(mz)) == 1) {
    if (verbose) {
      cat("Matching the m/z value...\n")
    }
    out <- "single"
  } else if (is.numeric(mz) && length(unique(mz)) == 2) {
    if (verbose) {
      cat("Matching the m/z interval...\n")
    }
    out <- "interval"
  } else if (is.numeric(mz) && length(unique(mz)) > 2) {
    if (verbose) {
      cat("Matching the", length(unique(mz)), "m/z values...\n")
    }
    out <- "array"
  } else {
    stop("Accepted values for mz are: 2 values vector for
         interval, array of values, or single mass value.")
  }

  return(out)
}

## .match.mz.single
.match.mz.single <- function(mz.value, mz.vector, mz.tol) {
  mz.vector <- sort(mz.vector)

  if (mz.value < mz.vector[1] || mz.value > mz.vector[length(mz.vector)]) {
    stop("M/Z outside the range of matched peaks m/z values.")
  }

  matched.idx <- c()
  for (i in 1:length(mz.vector))
  {
    if (.distance.ppm(mz.vector[i], mz.value) <= mz.tol) {
      matched.idx <- c(matched.idx, i)
    }
  }

  if (length(matched.idx) == 0) {
    stop("No m/z matched.")
  }

  if (length(matched.idx) > 1) {
    warning("Multiple matched m/z. Selecting the closest m/z value.\n")
    d <- abs(mz.vector[matched.idx] - mz.value)
    matched.idx <- matched.idx[which.min(d)]
  }

  cat("Matched m/z =", mz.vector[matched.idx], ".\n")

  return(matched.idx)
}

## .match.mz.interval
.match.mz.interval <- function(mz.values, mz.vector, mz.tol) {
  mz.values <- sort(mz.values)
  mz.vector <- sort(mz.vector)

  if (mz.tol <= 0) {
    stop("mz.tol must be positive")
  }

  if (mz.values[1] < mz.vector[1] || mz.values[2] > mz.vector[length(mz.vector)]) {
    stop("m/z range provided lays outside the measured m/z vector.")
  }

  # Check the closest masses within the mz.hw.ppm half window
  matched.lx <- c()
  matched.rx <- c()
  for (i in 1:length(mz.vector))
  {
    if (.distance.ppm(mz.vector[i], mz.values[1]) <= mz.tol) {
      matched.lx <- c(matched.lx, i)
    }
    if (.distance.ppm(mz.vector[i], mz.values[2]) <= mz.tol) {
      matched.rx <- c(matched.rx, i)
    }
  }
  if (length(matched.lx) == 0) {
    stop("No m/z matched the first element of the m/z interval.")
  }
  if (length(matched.rx) == 0) {
    stop("No m/z matched the last element of the m/z interval.")
  }
  # In case of multiple matched peaks, select the closest
  if (length(matched.lx) > 1) {
    warning("Multiple matched m/z. Selecting the closest m/z value.")
    d <- abs(mz.vector[matched.lx] - mz.values[1])
    matched.lx <- matched.lx[which.min(d)]
  }
  if (length(matched.rx) > 1) {
    warning("Multiple matched m/z. Selecting the closest m/z value.")
    d <- abs(mz.vector[matched.rx] - mz.values[2])
    matched.rx <- matched.rx[which.min(d)]
  }

  cat("Matched m/z endpoints:", mz.vector[matched.lx], ",", mz.vector[matched.rx], ".\n")
  cat("Num. selected peaks:", length(mz.vector[matched.lx:matched.rx]), ".\n")

  # Return the matched endpoints
  return(c(matched.lx, matched.rx))
}

## .match.mz.array
.match.mz.array <- function(mz.values, mz.vector, mz.tol) {
  mz.values <- sort(mz.values)
  mz.vector <- sort(mz.vector)

  # If some of the input mz.values are outside the range of mz.vector, stop.
  if (any(mz.values < mz.vector[1]) || any(mz.values > mz.vector[length(mz.vector)])) {
    stop("Some m/z outside the range of matched peaks m/z values.")
  }

  matched.idx <- vector(mode = "list", length = length(mz.values))
  for (i in 1:length(mz.values))
  {
    for (j in 1:length(mz.vector))
    {
      if (.distance.ppm(mz.values[i], mz.vector[j]) <= mz.tol) {
        matched.idx[[i]] <- c(matched.idx[[i]], j)
      }
    }
    rm(j)
  }
  rm(i)
  # Check whether some mz.values were matched with multiple values of mz.vector
  # or none of them
  l <- unlist(lapply(matched.idx, length))
  if (any(l == 0)) {
    stop("M/Z not matched: ", paste0(mz.values[l == 0], collapse = ", "), ".")
  }
  # These are the mz.values matched with 1 value of mz.vector
  fin.mz.idx <- unlist(matched.idx[l == 1])
  # If there are some matched with more than 1
  if (any(l > 1)) {
    warning("Multiple matched m/z. Selecting the closest m/z value.")
    multi_ <- which(l > 1)
    matched.idx <- matched.idx[multi_]
    # Now matched.idx contains only the mz.values matched with multiple elements
    # of mz.vector. For each of those elements, take the closest.
    for (i in 1:length(matched.idx))
    {
      # Calculate the distances between the mz.value[multi_[i]] and the matched
      # elements in mz.vector
      d <- numeric(length(matched.idx[[i]]))
      for (j in 1:length(matched.idx[i]))
      {
        d[j] <- abs(mz.vector[matched.idx[[i]][j]] - mz.values[multi_[i]])
      }
      rm(j)
      fin.mz.idx <- c(fin.mz.idx, matched.idx[[i]][which.min(d)])
    }
    rm(i)
  }

  cat("Matched m/z =", paste0(mz.vector[fin.mz.idx], collapse = ", "), ".\n")

  # Return the matched indices
  return(fin.mz.idx)
}
