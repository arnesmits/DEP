context("3 - Missing data imputation")

test_that("se2msn throws error without valid input", {
  expect_error(se2msn("test_se"))
})

test_that("se2msn returns a MSnSet object", {
  expect_is(se2msn(test_se), "MSnSet")
})

test_that("manual_impute throws error without valid input", {
  expect_error(manual_impute("test_vsn", 0.3, 1.8))
  expect_error(manual_impute(test_vsn, "0.3", 1.8))
  expect_error(manual_impute(test_vsn, 0.3, "1.8"))
})

test_that("manual_impute returns a MSnSet object", {
  expect_is(manual_impute(test_vsn), "SummarizedExperiment")
})

test_that("manual_impute returns an object without missing values", {
  result <- SummarizedExperiment::assay(manual_impute(test_vsn))
  expect_true(all(!is.na(result)))
})

test_that("imputation throws error without valid input", {
  expect_error(imputation("test_vsn", "MinProb"))
  expect_error(imputation(test_vsn, MinProb))
  expect_error(imputation(test_vsn, "FOO"))

  test_vsn_error <- test_vsn
  SummarizedExperiment::rowData(test_vsn_error) <- SummarizedExperiment::rowData(test_vsn_error)[,-(24:25)]
  expect_error(imputation(test_vsn_error, "MinProb"))
})

test_that("imputation returns a MSnSet object", {
  expect_is(imputation(test_vsn, "MinProb"), "SummarizedExperiment")
  expect_is(imputation(test_vsn, "man"), "SummarizedExperiment")
})

test_that("linear_model returns an object without missing values", {
  result <- SummarizedExperiment::assay(imputation(test_vsn, "MinProb"))
  expect_true(all(!is.na(result)))
})
