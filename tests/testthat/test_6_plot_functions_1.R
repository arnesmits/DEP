context("6 - Plot functions_1")

test_that("plot_single throws error without valid input", {
  expect_error(plot_single("test_sign", "USP19", "centered"))
  expect_error(plot_single(test_sign, USP19, "centered"))
  expect_error(plot_single(test_sign, "USP19", centered))
  expect_error(plot_single(test_sign, "USP19", "test"))
  expect_error(plot_single(test_sign, "USP15", "centered"))
  expect_error(plot_single(test_sign, "test", "centered"))

  test_sign_error1 <- test_sign
  SummarizedExperiment::colData(test_sign_error1) <- SummarizedExperiment::colData(test_sign_error1)[,-(3)]
  expect_error(plot_single(test_sign_error1, "USP19", "centered"))

  test_sign_error2 <- test_sign
  SummarizedExperiment::rowData(test_sign_error2) <- SummarizedExperiment::rowData(test_sign_error2)[,-(1)]
  expect_error(plot_single(test_sign_error2, "USP19", "centered"))

  test_sign_error3 <- test_sign
  SummarizedExperiment::rowData(test_sign_error3) <- SummarizedExperiment::rowData(test_sign_error3)[,-(26:31)]
  expect_error(plot_single(test_sign_error3, "USP19", "centered"))
})

test_that("plot_single returns a ggplot object", {
  expect_is(plot_single(test_sign, "USP19", "centered"), "ggplot")
  expect_is(plot_single(test_sign, "USP19", "contrast"), "ggplot")
})

test_that("plot_heatmap throws error without valid input", {
  expect_error(plot_heatmap("test_sign", "centered", k = 6, col_limit = 4, labelsize = 8))
  expect_error(plot_heatmap(test_sign, centered, k = 6, col_limit = 4, labelsize = 8))
  expect_error(plot_heatmap(test_sign, "centered", k = "6", col_limit = 4, labelsize = 8))
  expect_error(plot_heatmap(test_sign, "centered", k = 6, col_limit = "4", labelsize = 8))
  expect_error(plot_heatmap(test_sign, "centered", k = 6, col_limit = 4, labelsize = "8"))
  expect_error(plot_heatmap(test_sign, "test", k = 6, col_limit = 4, labelsize = 8))

  test_sign_error1 <- test_sign
  SummarizedExperiment::colData(test_sign_error1) <- SummarizedExperiment::colData(test_sign_error1)[,-(3)]
  expect_error(plot_heatmap(test_sign_error1, "centered", k = 6, col_limit = 4, labelsize = 8))

  test_sign_error2 <- test_sign
  SummarizedExperiment::rowData(test_sign_error2) <- SummarizedExperiment::rowData(test_sign_error2)[,-(35)]
  expect_error(plot_heatmap(test_sign_error2, "centered", k = 6, col_limit = 4, labelsize = 8))

  test_sign_error3 <- test_sign
  SummarizedExperiment::rowData(test_sign_error3) <- SummarizedExperiment::rowData(test_sign_error3)[,-(26:28)]
  expect_error(plot_heatmap(test_sign_error3, "centered", k = 6, col_limit = 4, labelsize = 8))
})

test_that("plot_heatmap returns a HeatmapList object", {
  expect_is(plot_heatmap(test_sign, "centered", k = 6, col_limit = 4, labelsize = 3), "HeatmapList")
  expect_is(plot_heatmap(test_sign, "contrast", k = 6, col_limit = 8, labelsize = 3), "HeatmapList")
  expect_is(plot_heatmap(test_ibaq_sign, "contrast", k = 1, col_limit = 8, labelsize = 3), "HeatmapList")
})

test_that("plot_volcano throws error without valid input", {
  expect_error(plot_volcano("test_sign", "Ubi6_vs_Ctrl", labelsize = 5, add_names = TRUE))
  expect_error(plot_volcano(test_sign, Ubi6_vs_Ctrl, labelsize = 5, add_names = TRUE))
  expect_error(plot_volcano(test_sign, "Ubi6_vs_Ctrl", labelsize = 5, add_names = "TRUE"))
  expect_error(plot_volcano(test_sign, "Ubi6_vs_Ctrl", labelsize = "5", add_names = TRUE))
  expect_error(plot_volcano(test_sign, "Ubi5_vs_Ctrl", labelsize = 5, add_names = TRUE))

  test_sign_error1 <- test_sign
  SummarizedExperiment::rowData(test_sign_error1) <- SummarizedExperiment::rowData(test_sign_error1)[,-(1)]
  expect_error(plot_volcano(test_sign_error1, "Ubi6_vs_Ctrl", labelsize = 5, add_names = TRUE))

  test_sign_error2 <- test_sign
  SummarizedExperiment::rowData(test_sign_error2) <- SummarizedExperiment::rowData(test_sign_error2)[,-(26:31)]
  expect_error(plot_volcano(test_sign_error2, "Ubi6_vs_Ctrl", labelsize = 5, add_names = TRUE))

  test_sign_error3 <- test_sign
  SummarizedExperiment::rowData(test_sign_error3) <- SummarizedExperiment::rowData(test_sign_error3)[,-(32:34)]
  expect_error(plot_volcano(test_sign_error3, "Ubi6_vs_Ctrl", labelsize = 5, add_names = TRUE))
})

test_that("plot_volcano returns a ggplot object", {
  expect_is(plot_volcano(test_sign, "Ubi6_vs_Ctrl", labelsize = 5, add_names = TRUE), "ggplot")
  expect_is(plot_volcano(test_sign, "Ubi6_vs_Ctrl", add_names = FALSE), "ggplot")
})
