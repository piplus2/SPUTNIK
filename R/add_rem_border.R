# Add a border to an image
addBorderImage <- function(imMat, border = 2) {
  if (!is.matrix(imMat)) {
    stop("input must be a numeric matrix.")
  }
  if (any(!is.finite(c(imMat)))) {
    stop("input must be finite.")
  }
  if (!is.finite(border)) {
    stop("border must be finite.")
  }
  if (border == 0) {
    return(imMat)
  }
  if (border < 0) {
    stop("border must be positive.")
  }

  imMat <- rbind(
    matrix(0, border, ncol(imMat)),
    imMat,
    matrix(0, border, ncol(imMat))
  )
  imMat <- cbind(
    matrix(0, nrow(imMat), border),
    imMat,
    matrix(0, nrow(imMat), border)
  )
  return(imMat)
}

# Remove the border from an image
remBorderImage <- function(imMat, border = 2) {
  if (!is.matrix(imMat)) {
    stop("input must be a numeric matrix.")
  }
  if (!all(is.finite(c(imMat)))) {
    stop("input must be finite.")
  }
  if (!is.finite(border)) {
    stop("border must be finite.")
  }
  if (border == 0) {
    return(imMat)
  }
  if (border < 0) {
    stop("border must be positive.")
  }

  imMat <- imMat[
    seq(border + 1, nrow(imMat) - border),
    seq(border + 1, ncol(imMat) - border)
  ]
  return(imMat)
}
