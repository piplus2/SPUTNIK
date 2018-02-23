## Load package
library("SPUTNIK")

## Image
im <- matrix(rnorm(100, 10, 10))

## Spatial chaos
sc <- spatial.chaos(im, levels = 30, morph = TRUE)

## Gini index
gi <- gini.index(im, levels = 16)

## Scatter ratio
sr <- scatter.ratio(im)
