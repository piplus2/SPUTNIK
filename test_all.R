library(SPUTNIK)
library(testthat)

source('tests/test_constructors.R')

source('tests/test_ext_data.R')

source('tests/test_global_ref_filter.R')

source('tests/test_pixel_count_filter.R')

source('tests/test_csr_filter.R')

test_results <- test_dir('tests/', reporter = 'summary')
print(test_results)