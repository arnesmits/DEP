context("9 - Plot functions explore")

test_that("plot_pca throws error without valid input", {
  expect_error(plot_pca("test_sign", x = 1, y = 2, n = 100))
  expect_error(plot_pca(test_sign, x = "1", y = 2, n = 100))
  expect_error(plot_pca(test_sign, x = 1, y = "2", n = 100))
  expect_error(plot_pca(test_sign, x = 1, y = 2, n = "100"))
  expect_error(plot_pca(test_sign, x = 20, y = 2, n = 100))
  expect_error(plot_pca(test_sign, x = 1, y = 2, n = 500))
})

test_that("plot_pca returns a ggplot object", {
  expect_is(plot_pca(test_sign, x = 1, y = 2, n = 100), "ggplot")
})

test_that("plot_corr throws error without valid input", {
  expect_error(plot_corr("test_sign", TRUE, -1, 1, "PRGn", FALSE))
  expect_error(plot_corr(test_sign, "TRUE", -1, 1, "PRGn", FALSE))
  expect_error(plot_corr(test_sign, TRUE, "-1", 1, "PRGn", FALSE))
  expect_error(plot_corr(test_sign, TRUE, -1, "1", "PRGn", FALSE))
  expect_error(plot_corr(test_sign, TRUE, -1, 1, PRGn, FALSE))
  expect_error(plot_corr(test_sign, TRUE, -1, 1, "PRGn", "FALSE"))
  expect_error(plot_corr(test_sign, TRUE, -1, 1, "Accent", FALSE))
  expect_error(plot_corr(test_sign, TRUE, -1, 2, "Accent", FALSE))

  test_sign_error <- test_sign
  SummarizedExperiment::rowData(test_sign_error) <- SummarizedExperiment::rowData(test_sign_error)[,-(35)]
  expect_error(plot_corr(test_sign_error, TRUE, -1, 1, "PRGn", FALSE))
})

test_that("plot_corr_matrix returns a ComplexHeatmap object", {
  expect_is(plot_corr(test_sign, TRUE, -1, 1, "PRGn", FALSE), "HeatmapList")
})
