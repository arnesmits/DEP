context("4 - Differential testing")

test_that("linear_model throws error without valid input", {
  expect_error(linear_model("test_impute", "Ctrl", "all"))
  expect_error(linear_model(test_impute, Ctrl, "all"))
  expect_error(linear_model(test_impute, "Ctrl", all))
  expect_error(linear_model(test_impute, "Ctrl", "ALL"))
  expect_error(linear_model(test_impute, "control", "all"))

  test_impute_error <- test_impute
  SummarizedExperiment::colData(test_impute_error) <- SummarizedExperiment::colData(test_impute_error)[,-(3)]
  expect_error(linear_model(test_impute_error, "Ctrl", "all"))

  test_impute_error2 <- test_impute
  SummarizedExperiment::rowData(test_impute_error2) <- SummarizedExperiment::rowData(test_impute_error2)[,-(24:25)]
  expect_error(linear_model(test_impute_error2, "Ctrl", "all"))
})

test_that("linear_model returns a SummarizedExperiment object", {
  expect_is(linear_model(test_impute, "Ctrl", "all"), "SummarizedExperiment")
})

test_that("linear_model returns an object with diff and p.adj columns", {
  result <- SummarizedExperiment::rowData(linear_model(test_impute, "Ctrl", "control"))
  expect_equal(grep("_diff$", colnames(result)), 26:28)
  expect_equal(grep("_p.adj$", colnames(result)), 29:31)
})

test_that("cutoffs throws error without valid input", {
  expect_error(cutoffs("test_lm", 0.05, 1))
  expect_error(cutoffs(test_lm, "0.05", 1))
  expect_error(cutoffs(test_lm, 0.05, "1"))
  expect_error(cutoffs(test_impute, 0.05, 1))

  test_lm_error <- test_lm
  SummarizedExperiment::rowData(test_lm_error) <- SummarizedExperiment::rowData(test_lm_error)[,-(24:25)]
  expect_error(cutoffs(test_lm_error, 0.05, 1))
})

test_that("cutoffs returns a SummarizedExperiment object", {
  expect_is(cutoffs(test_lm, 0.05, 1), "SummarizedExperiment")
})

test_that("cutoffs returns an object with significance columns", {
  result <- SummarizedExperiment::rowData(cutoffs(test_lm, 0.05, 1))
  expect_equal(grep("^sign$", colnames(result)), 35)
  expect_equal(grep("_sign$", colnames(result)), 32:34)
  expect_equal(nrow(result[result$sign == "+",]), 45)
})

test_that("results throws error without valid input", {
  expect_error(results("test_sign"))
  expect_error(results(test_impute))

  test_sign_error <- test_sign
  SummarizedExperiment::rowData(test_sign_error) <- SummarizedExperiment::rowData(test_sign_error)[,-(24:25)]
  expect_error(results(test_sign_error))
})

test_that("results returns a data.frame", {
  expect_is(results(test_sign), "data.frame")
})

test_that("results returns a data.frame with the expected column", {
  result <- results(test_sign)
  expect_equal(nrow(result), 203)
  expect_equal(ncol(result), 16)
  expect_equal(grep("_centered$", colnames(result)), 13:16)
  expect_equal(grep("_ratio$", colnames(result)), 10:12)
  expect_equal(nrow(result[result$sign == "+",]), 45)
})
