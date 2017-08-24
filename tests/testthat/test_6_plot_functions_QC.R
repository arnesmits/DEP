context("7 - Plot functions QC")

test_that("meanSdPlot returns a list object", {
  expect_is(meanSdPlot(test_vsn), "list")
})

test_that("plot_normalization throws error without valid input", {
  expect_error(plot_normalization("test_filter", test_vsn))
  expect_error(plot_normalization(test_filter, "test_vsn"))

  test_filter_error <- test_filter
  SummarizedExperiment::colData(test_filter_error) <- SummarizedExperiment::colData(test_filter_error)[,-(3)]
  expect_error(plot_normalization(test_filter_error, test_vsn))

  test_vsn_error <- test_vsn
  SummarizedExperiment::colData(test_vsn_error) <- SummarizedExperiment::colData(test_vsn_error)[,-(3)]
  expect_error(plot_normalization(test_filter, test_vsn_error))
})

test_that("plot_normalization returns a ggplot object", {
  expect_is(plot_normalization(test_filter, test_vsn), "ggplot")
})

test_that("plot_imputation throws error without valid input", {
  expect_error(plot_imputation("test_vsn", test_impute))
  expect_error(plot_imputation(test_vsn, "test_impute"))

  test_vsn_error <- test_vsn
  SummarizedExperiment::colData(test_vsn_error) <- SummarizedExperiment::colData(test_vsn_error)[,-(3)]
  expect_error(plot_imputation(test_vsn_error, test_impute))

  test_impute_error <- test_impute
  SummarizedExperiment::colData(test_impute_error) <- SummarizedExperiment::colData(test_impute_error)[,-(3)]
  expect_error(plot_imputation(test_vsn, test_impute_error))
})

test_that("plot_imputation returns a ggplot object", {
  expect_is(plot_imputation(test_vsn, test_impute), "ggplot")
})

test_that("plot_detect throws error without valid input", {
  expect_error(plot_detect("test_filter"))
})

test_that("plot_detect returns a grob object", {
  expect_is(plot_detect(test_filter), "grob")
})

test_that("plot_missval throws error without valid input", {
  expect_error(plot_missval("test_filter"))

  test_filter_noNA <-
    test_filter[!apply(assay(test_filter), 1, function(x) any(is.na(x))),]
  expect_error(plot_missval(test_filter_noNA))
})

test_that("plot_missval returns a HeatmapList object", {
  expect_is(plot_missval(test_filter), "HeatmapList")
})
