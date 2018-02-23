## .scaleAllVariables
.scale.all <- function(msi.data)
{
  .stopIfNotValidMSIDataset(msi.data)

  x <- msi.data@matrix
  for (i in 1:ncol(x))
  {
    x[, i] <- x[, i] / max(x[, i])
  }
  msi.data@matrix <- x

  msi.data
}
