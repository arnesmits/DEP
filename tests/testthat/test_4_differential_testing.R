context("4 - Differential testing")

test_that("test_diff throws error without valid input", {
  expect_error(test_diff("test_impute", "Ctrl", "all"))
  expect_error(test_diff(test_impute, Ctrl, "all"))
  expect_error(test_diff(test_impute, "Ctrl", all))
  expect_error(test_diff(test_impute, "Ctrl", "ALL"))
  expect_error(test_diff(test_impute, "control", "all"))
  expect_error(test_diff(test_impute, "Ctrl", "manual"))
  expect_error(test_diff(test_impute, "Ctrl", "manual", Ubi4_vs_Ctrl))
  expect_error(test_diff(test_impute, "Ctrl", "manual", "test"))
  expect_error(test_diff(test_impute, "Ctrl", "control", incl_repl = "TRUE"))

  test_impute_error <- test_impute
  SummarizedExperiment::colData(test_impute_error) <- SummarizedExperiment::colData(test_impute_error)[,-(3)]
  expect_error(test_diff(test_impute_error, "Ctrl", "all"))

  test_impute_error2 <- test_impute
  SummarizedExperiment::rowData(test_impute_error2) <- SummarizedExperiment::rowData(test_impute_error2)[,-(24:25)]
  expect_error(test_diff(test_impute_error2, "Ctrl", "all"))
})

test_that("test_diff returns a SummarizedExperiment object", {
  expect_is(test_diff(test_impute, "Ctrl", "all"), "SummarizedExperiment")
  expect_is(test_diff(test_impute, "Ctrl", "manual", "Ubi4_vs_Ctrl"), "SummarizedExperiment")
  expect_is(test_diff(test_impute, "Ctrl", "all", incl_repl = TRUE), "SummarizedExperiment")
})

test_that("test_diff returns an object with diff and p.adj columns", {
  result <- SummarizedExperiment::rowData(test_diff(test_impute, "Ctrl", "control"))
  expect_equal(grep("_diff$", colnames(result)), 26:28)
  expect_equal(grep("_p.adj$", colnames(result)), 29:31)
})

test_that("add_rejections throws error without valid input", {
  expect_error(add_rejections("test_lm", 0.05, 1))
  expect_error(add_rejections(test_lm, "0.05", 1))
  expect_error(add_rejections(test_lm, 0.05, "1"))
  expect_error(add_rejections(test_impute, 0.05, 1))

  test_lm_error <- test_lm
  SummarizedExperiment::rowData(test_lm_error) <- SummarizedExperiment::rowData(test_lm_error)[,-(24:25)]
  expect_error(add_rejections(test_lm_error, 0.05, 1))
})

test_that("add_rejections returns a SummarizedExperiment object", {
  expect_is(add_rejections(test_lm, 0.05, 1), "SummarizedExperiment")
})

test_that("add_rejections returns an object with significance columns", {
  result <- SummarizedExperiment::rowData(add_rejections(test_lm, 0.05, 1))
  expect_equal(grep("^significant$", colnames(result)), 35)
  expect_equal(grep("_significant$", colnames(result)), 32:34)
  expect_equal(nrow(result[result$significant,]), 45)
})

test_that("get_results throws error without valid input", {
  expect_error(get_results("test_sign"))
  expect_error(get_results(test_impute))

  test_sign_error <- test_sign
  SummarizedExperiment::rowData(test_sign_error) <- SummarizedExperiment::rowData(test_sign_error)[,-(24:25)]
  expect_error(get_results(test_sign_error))
})

test_that("get_results returns a data.frame", {
  expect_is(get_results(test_sign), "data.frame")
})

test_that("get_results returns a data.frame with the expected column", {
  result <- get_results(test_sign)
  expect_equal(nrow(result), 203)
  expect_equal(ncol(result), 16)
  expect_equal(grep("_centered$", colnames(result)), 13:16)
  expect_equal(grep("_ratio$", colnames(result)), 10:12)
  expect_equal(nrow(result[result$significant,]), 45)
})

test_that("get_df_wide throws error without valid input", {
  expect_error(get_df_wide("test_sign"))

  test_sign_error <- test_sign
  SummarizedExperiment::rowData(test_sign_error) <- SummarizedExperiment::rowData(test_sign_error)[,-(1)]
  expect_error(get_df_wide(test_sign_error))
})

test_that("get_df_wide returns a data.frame", {
  expect_is(get_df_wide(test_sign), "data.frame")
})

test_that("get_df_long throws error without valid input", {
  expect_error(get_df_long("test_sign"))

  test_sign_error <- test_sign
  SummarizedExperiment::rowData(test_sign_error) <- SummarizedExperiment::rowData(test_sign_error)[,-(1)]
  expect_error(get_df_long(test_sign_error))
})

test_that("get_df_long returns a data.frame", {
  expect_is(get_df_long(test_sign), "data.frame")
})
