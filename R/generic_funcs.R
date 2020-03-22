## .scaleAllVariables
.scale.all <- function(msi.data) {
  .stopIfNotValidMSIDataset(msi.data)

  x <- msi.data@matrix
  for (i in 1:ncol(x))
  {
    x[, i] <- x[, i] / max(x[, i])
  }
  msi.data@matrix <- x

  return(msi.data)
}

## Statistical mode
.mode <- function(x) {
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

## removeConstPeaks 
.remove.const.peaks = function(object) {
  .stopIfNotValidMSIDataset(object)
  const.filter <- apply(object@matrix, 2, var) == 0
  if (sum(const.filter) != 0) {
    cat(paste0("Found ", sum(const.filter), " constant peaks.\n"))
    object@matrix <- object@matrix[, !const.filter]
    object@mz <- object@mz[!const.filter]
  }
  return(object)
}
