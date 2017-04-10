context("7 - Plot functions2")

test_that("plot_norm throws error without valid input", {
  expect_error(plot_norm("test_filter", test_vsn))
  expect_error(plot_norm(test_filter, "test_vsn"))

  test_filter_error <- test_filter
  SummarizedExperiment::colData(test_filter_error) <- SummarizedExperiment::colData(test_filter_error)[,-(3)]
  expect_error(plot_norm(test_filter_error, test_vsn))

  test_vsn_error <- test_vsn
  SummarizedExperiment::colData(test_vsn_error) <- SummarizedExperiment::colData(test_vsn_error)[,-(3)]
  expect_error(plot_norm(test_filter, test_vsn_error))
})

test_that("plot_norm returns a ggplot object", {
  expect_is(plot_norm(test_filter, test_vsn), "ggplot")
})

test_that("plot_imp throws error without valid input", {
  expect_error(plot_imp("test_vsn", test_impute))
  expect_error(plot_imp(test_vsn, "test_impute"))

  test_vsn_error <- test_vsn
  SummarizedExperiment::colData(test_vsn_error) <- SummarizedExperiment::colData(test_vsn_error)[,-(3)]
  expect_error(plot_imp(test_vsn_error, test_impute))

  test_impute_error <- test_impute
  SummarizedExperiment::colData(test_impute_error) <- SummarizedExperiment::colData(test_impute_error)[,-(3)]
  expect_error(plot_imp(test_vsn, test_impute_error))
})

test_that("plot_imp returns a ggplot object", {
  expect_is(plot_imp(test_vsn, test_impute), "ggplot")
})

test_that("plot_detect throws error without valid input", {
  expect_error(plot_detect("test_filter"))

  test_filter_error <- test_filter
  SummarizedExperiment::colData(test_filter_error) <- SummarizedExperiment::colData(test_filter_error)[,-(3)]
  expect_error(plot_detect(test_filter_error))
})

test_that("plot_detect returns a grob object", {
  expect_is(plot_detect(test_filter), "grob")
})

test_that("plot_missval throws error without valid input", {
  expect_error(plot_missval("test_filter"))
})

test_that("plot_missval returns a HeatmapList object", {
  expect_is(plot_missval(test_filter), "HeatmapList")
})

test_that("plot_numbers throws error without valid input", {
  expect_error(plot_numbers("test_filter"))
})

test_that("plot_numbers returns a ggplot object", {
  expect_is(plot_numbers(test_filter), "ggplot")
})

test_that("plot_coverage throws error without valid input", {
  expect_error(plot_coverage("test_filter"))
})

test_that("plot_coverage returns a ggplot object", {
  expect_is(plot_coverage(test_filter), "ggplot")
})

test_that("plot_frequency throws error without valid input", {
  expect_error(plot_frequency("test_filter"))
})

test_that("plot_frequency returns a ggplot object", {
  expect_is(plot_frequency(test_filter), "ggplot")
})
