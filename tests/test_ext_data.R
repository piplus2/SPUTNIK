# Test external data

library(testthat)
library(SPUTNIK)

# MALDI-MSI example

test_that("bladder MALDI-MSI", {
  x <- bladderMALDIRompp2010(verbose = TRUE)
  expect_equal(length(attr(x, "mass")), 1175)
  expect_equal(attr(x, "size"), c(260, 134))
  expect_equal(dim(x), c(34840, 1175))
})

test_that("ovarian DESI-MSI", {
  x <- ovarianDESIDoria2016(verbose = TRUE)
  expect_equal(length(attr(x, "mass")), 1375)
  expect_equal(attr(x, "size"), c(73, 118))
  expect_equal(dim(x), c(8614, 1375))
})
