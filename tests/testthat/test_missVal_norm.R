context("2 - Missing values filtering, data normalization and missing values imputation")

test_that("filter_missval throws error without valid input", {
  expect_error(filter_missval("test_se", 0))
  expect_error(filter_missval(test_se, "0"))
})

test_that("filter_missval returns a SummarizedExperiment", {
  expect_is(filter_missval(test_se, 0), "SummarizedExperiment")
})

test_that("filter_missval returns correct number of rows", {
  expect_equal(filter_missval(test_se, 0) %>% nrow(), 203)
  expect_equal(filter_missval(test_se, 1) %>% nrow(), 243)
  expect_equal(filter_missval(test_se, 2) %>% nrow(), 293)
})
