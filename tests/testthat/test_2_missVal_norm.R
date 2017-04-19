context("2 - Missing values filtering and data normalization")

test_that("filter_missval throws error without valid input", {
  expect_error(filter_missval("test_se", 0))
  expect_error(filter_missval(test_se, "0"))
  expect_error(filter_missval(test_se, 4))
  expect_error(filter_missval(test_se, -1))

  test_se_error <- test_se
  SummarizedExperiment::rowData(test_se_error) <- SummarizedExperiment::rowData(test_se_error)[,-(24:25)]
  expect_error(filter_missval(test_se_error, 0))

  test_se_error2 <- test_se
  SummarizedExperiment::colData(test_se_error2) <- SummarizedExperiment::colData(test_se_error2)[,-(3)]
  expect_error(filter_missval(test_se_error2, 0))
})

test_that("filter_missval returns a SummarizedExperiment", {
  expect_is(filter_missval(test_se, 0), "SummarizedExperiment")
})

test_that("filter_missval returns correct number of rows", {
  expect_equal(filter_missval(test_se, 0) %>% nrow(), 203)
  expect_equal(filter_missval(test_se, 1) %>% nrow(), 243)
  expect_equal(filter_missval(test_se, 2) %>% nrow(), 393)
})

test_that("normalize throws error without valid input", {
  expect_error(normalize("test_filter"))
})

test_that("normalize returns a SummarizedExperiment", {
  expect_is(normalize(test_filter), "SummarizedExperiment")
})
