context("5 - iBAQ functions")

test_that("merge_ibaq throws error without valid input", {
  expect_error(merge_ibaq("test_ibaq", test_ibaq_pep))
  expect_error(merge_ibaq(test_ibaq, "test_ibaq_pep"))
  expect_error(merge_ibaq(test_ibaq[,-(30)], test_ibaq_pep))
  expect_error(merge_ibaq(test_ibaq[,-(31)], test_ibaq_pep))
  expect_error(merge_ibaq(test_ibaq[,-(15:20)], test_ibaq_pep))
  expect_error(merge_ibaq(test_ibaq, test_ibaq_pep[,-(14)]))
  expect_error(merge_ibaq(test_ibaq, test_ibaq_pep[,-(6)]))
})

test_that("merge_ibaq returns a data.frame", {
  expect_is(merge_ibaq(test_ibaq, test_ibaq_pep), "data.frame")
})

test_that("merge_ibaq returns an object with the rigth dimensions and columns", {
  result <- merge_ibaq(test_ibaq, test_ibaq_pep)
  expect_equal(grep("iBAQ", colnames(result)), 4:9)
  expect_equal(dim(result), c(359,10))
})

test_that("stoichiometry throws error without valid input", {
  expect_error(stoichiometry("test_ibaq_sign", test_ibaq_ibaq, "GFP_vs_WT", "Rbbp4", 1))
  expect_error(stoichiometry(test_ibaq_sign, "test_ibaq_ibaq", "GFP_vs_WT", "Rbbp4", 1))
  expect_error(stoichiometry(test_ibaq_sign, test_ibaq_ibaq, GFP_vs_WT, "Rbbp4", 1))
  expect_error(stoichiometry(test_ibaq_sign, test_ibaq_ibaq, "GFP_vs_WT", Rbbp4, 1))
  expect_error(stoichiometry(test_ibaq_sign, test_ibaq_ibaq, "GFP_vs_WT", "Rbbp4", "1"))
  expect_error(stoichiometry(test_ibaq_sign, test_ibaq_ibaq[,-(2:3)], "GFP_vs_WT", "Rbbp4", 1))
  expect_error(stoichiometry(test_ibaq_sign, test_ibaq_ibaq[,-(4:9)], "GFP_vs_WT", "Rbbp4", 1))

  test_ibaq_sign_error1 <- test_ibaq_sign
  SummarizedExperiment::rowData(test_ibaq_sign_error1) <- SummarizedExperiment::rowData(test_ibaq_sign_error1)[,-(1)]
  expect_error(stoichiometry(test_ibaq_sign_error1, test_ibaq_ibaq, "GFP_vs_WT", "Rbbp4", 1))

  test_ibaq_sign_error2 <- test_ibaq_sign
  SummarizedExperiment::rowData(test_ibaq_sign_error2) <- SummarizedExperiment::rowData(test_ibaq_sign_error2)[,-c(27,28,30)]
  expect_error(stoichiometry(test_ibaq_sign_error2, test_ibaq_ibaq, "GFP_vs_WT", "Rbbp4", 1))

})

test_that("stoichiometry returns a data.frame", {
  expect_is(stoichiometry(test_ibaq_sign, test_ibaq_ibaq, "GFP_vs_WT", "Rbbp4", 1), "data.frame")
})

test_that("stoichiometry returns an object with the rigth dimensions and columns", {
  result <- stoichiometry(test_ibaq_sign, test_ibaq_ibaq, "GFP_vs_WT", "Rbbp4", 1)
  expect_equal(result$stoichiometry[result$name == "Rbbp4"], 1)
  expect_equal(dim(result), c(5,4))
})

test_that("plot_stoi throws error without valid input", {
  expect_error(plot_stoi("test_stoi", 0.001, NULL))
  expect_error(plot_stoi(test_stoi, "0.001", NULL))
  expect_error(plot_stoi(test_stoi, 0.001, "0.05"))
  expect_error(plot_stoi(test_stoi[,-(1)], 0.001, NULL))
  expect_error(plot_stoi(test_stoi[,-(2)], 0.001, NULL))
  expect_error(plot_stoi(test_stoi[,-(3)], 0.001, NULL))
  expect_error(plot_stoi(test_stoi[,-(4)], 0.001, NULL))
})

test_that("plot_stoi returns a ggplot object", {
  expect_is(plot_stoi(test_stoi, 0.001), "ggplot")
  expect_is(plot_stoi(test_stoi, 0.001, 0.5), "ggplot")
})

test_that("plot_ibaq throws error without valid input", {
  expect_error(plot_ibaq("test_ibaq_sign", "GFP_vs_WT", 3))
  expect_error(plot_ibaq(test_ibaq_sign, GFP_vs_WT, 3))
  expect_error(plot_ibaq(test_ibaq_sign, "GFP_vs_WT", "3"))
  expect_error(plot_ibaq(test_ibaq_sign, "test", 3))

  test_ibaq_sign_error1 <- test_ibaq_sign
  SummarizedExperiment::rowData(test_ibaq_sign_error1) <- SummarizedExperiment::rowData(test_ibaq_sign_error1)[,-(1)]
  expect_error(plot_ibaq(test_ibaq_sign_error1, "GFP_vs_WT", 3))

  test_ibaq_sign_error2 <- test_ibaq_sign
  SummarizedExperiment::rowData(test_ibaq_sign_error2) <- SummarizedExperiment::rowData(test_ibaq_sign_error2)[,-c(27,28,30)]
  expect_error(plot_ibaq(test_ibaq_sign_error2, "GFP_vs_WT", 3))

  test_ibaq_sign_error3 <- test_ibaq_sign
  SummarizedExperiment::rowData(test_ibaq_sign_error3) <- SummarizedExperiment::rowData(test_ibaq_sign_error3)[,-(16:21)]
  expect_error(plot_ibaq(test_ibaq_sign_error3, "GFP_vs_WT", 3))

})

test_that("plot_ibaq returns a ggplot object", {
  expect_is(plot_ibaq(test_ibaq_sign, "GFP_vs_WT", 3), "ggplot")
})
