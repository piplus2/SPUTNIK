## Load package
library("SPUTNIK")

## Create the msi.dataset-class object
sz <- c(40, 40)
x <- matrix(rnorm(sz[1] * sz[2] * 20) * 1000, sz[1] * sz[2], 20)
x[x < 0] <- 0 # MS data is positive
mz <- sort(sample(100, ncol(x)))
msiX <- msiDataset(x, mz, sz[1], sz[2])

## Normalize and log-transform
msiX <- normIntensity(msiX, "median")
msiX <- varTransform(msiX, "log")

## Create the msi.dataset-class object
sz <- c(40, 40)
x <- matrix(rnorm(sz[1] * sz[2] * 20) * 1000, sz[1] * sz[2], 20)
x[x < 0] <- 0 # MS data is positive
mz <- sort(sample(100, ncol(x)))
msiX <- msiDataset(x, mz, sz[1], sz[2])

## Normalize using PQN
msiX <- normIntensity(msiX, "PQN")
