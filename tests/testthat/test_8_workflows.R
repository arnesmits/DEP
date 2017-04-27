context("8 - Workflows")

test_that("LFQ throws error without valid input", {
  expect_error(LFQ("test_data", UbiLength_ExpDesign, "MinProb", "Ctrl", "control", filter = c("Reverse", "Potential.contaminant"),
                   name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1))
  expect_error(LFQ(test_data, "UbiLength_ExpDesign", "MinProb", "Ctrl", "control", filter = c("Reverse", "Potential.contaminant"),
                   name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1))
  expect_error(LFQ(test_data, UbiLength_ExpDesign, MinProb, "Ctrl", "control", filter = c("Reverse", "Potential.contaminant"),
                   name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1))
  expect_error(LFQ(test_data, UbiLength_ExpDesign, "MinProb", Ctrl, "control", filter = c("Reverse", "Potential.contaminant"),
                   name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1))
  expect_error(LFQ(test_data, UbiLength_ExpDesign, "MinProb", "Ctrl", control, filter = c("Reverse", "Potential.contaminant"),
                   name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1))
  expect_error(LFQ(test_data, UbiLength_ExpDesign, "MinProb", "Ctrl", "control", filter = test,
                   name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1))
  expect_error(LFQ(test_data, UbiLength_ExpDesign, "MinProb", "Ctrl", "control", filter = c("Reverse", "Potential.contaminant"),
                   name = Gene.names, ids = "Protein.IDs", alpha = 0.05, lfc = 1))
  expect_error(LFQ(test_data, UbiLength_ExpDesign, "MinProb", "Ctrl", "control", filter = c("Reverse", "Potential.contaminant"),
                   name = "Gene.names", ids = Protein.IDs, alpha = 0.05, lfc = 1))
  expect_error(LFQ(test_data, UbiLength_ExpDesign, "MinProb", "Ctrl", "control", filter = c("Reverse", "Potential.contaminant"),
                   name = "Gene.names", ids = "Protein.IDs", alpha = "0.05", lfc = 1))
  expect_error(LFQ(test_data, UbiLength_ExpDesign, "MinProb", "Ctrl", "control", filter = c("Reverse", "Potential.contaminant"),
                   name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = "1"))

  expect_error(LFQ(test_data, UbiLength_ExpDesign, "MinProb", "Ctrl", "control", filter = "test",
                   name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1))
  expect_error(LFQ(test_data, UbiLength_ExpDesign, "Test", "Ctrl", "control", filter = c("Reverse", "Potential.contaminant"),
                   name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1))
  expect_error(LFQ(test_data, UbiLength_ExpDesign, "MinProb", "Ctrl", "control", filter = c("Reverse", "Potential.contaminant"),
                   name = "Gene.Names", ids = "Protein.IDs", alpha = 0.05, lfc = 1))
  expect_error(LFQ(test_data, UbiLength_ExpDesign, "MinProb", "Ctrl", "control", filter = c("Reverse", "Potential.contaminant"),
                   name = "Gene.names", ids = "Protein.ids", alpha = 0.05, lfc = 1))
  expect_error(LFQ(test_data, UbiLength_ExpDesign[,-(2)], "MinProb", "Ctrl", "control", filter = c("Reverse", "Potential.contaminant"),
                   name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1))
  expect_error(LFQ(test_data, UbiLength_ExpDesign, "MinProb", "ctrl", "control", filter = c("Reverse", "Potential.contaminant"),
                   name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1))
  expect_error(LFQ(test_data, UbiLength_ExpDesign, "MinProb", "Ctrl", "test", filter = c("Reverse", "Potential.contaminant"),
                   name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1))
})

test_that("LFQ returns a list object", {
  expect_is(LFQ(test_data, UbiLength_ExpDesign, "MinProb", "Ctrl", "control", filter = c("Reverse", "Potential.contaminant"),
                name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1), "list")
  expect_is(LFQ(test_data, UbiLength_ExpDesign, "MinProb", "Ctrl", "control", filter = c("Reverse", "Potential.contaminant"),
                name = "Gene.names", ids = "Protein.IDs", alpha = 0.05L, lfc = 1L), "list")
  expect_is(LFQ(test_data, UbiLength_ExpDesign, "MinProb", "Ctrl", "control", filter = "Reverse",
                name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1), "list")
})

test_that("LFQ returns a list object with all expected objects", {
  expect_equal(names(LFQ(test_data, UbiLength_ExpDesign, "MinProb", "Ctrl", "control", filter = c("Reverse", "Potential.contaminant"),
                name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1)), c("data","se","filt","norm","imputed","diff","signif","results","param"))
  expect_equal(names(LFQ(tibble::as_tibble(test_data), UbiLength_ExpDesign, "MinProb", "Ctrl", "control", filter = c("Reverse", "Potential.contaminant"),
                         name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1)), c("data","se","filt","norm","imputed","diff","signif","results","param"))
  expect_equal(names(LFQ(test_data, tibble::as_tibble(UbiLength_ExpDesign), "MinProb", "Ctrl", "control", filter = c("Reverse", "Potential.contaminant"),
                         name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1)), c("data","se","filt","norm","imputed","diff","signif","results","param"))
})

test_that("report throws error without valid input", {
  res <- LFQ(test_data, UbiLength_ExpDesign, "MinProb", "Ctrl", "control", filter = c("Reverse", "Potential.contaminant"),
             name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1)
  expect_error(report("res"))
  expect_error(report(res[-1]))
})

test_that("report output ", {
  res <- LFQ(test_data, UbiLength_ExpDesign, "MinProb", "Ctrl", "control", filter = c("Reverse", "Potential.contaminant"),
             name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1)
  expect_output(report(res), "Report and results.txt file saved in:")
})

test_that("iBAQ throws error without valid input", {
  result <- LFQ(test_ibaq[,-(31:32)], GFPip_ExpDesign, "MinProb", "WT", "control", filter = c("Reverse", "Contaminant"),
                name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1)
  expect_error(iBAQ("result", test_ibaq_pep, "GFP_vs_WT", "Rbbp4"))
  expect_error(iBAQ(result, "test_ibaq_pep", "GFP_vs_WT", "Rbbp4"))
  expect_error(iBAQ(result, test_ibaq_pep, GFP_vs_WT, "Rbbp4"))
  expect_error(iBAQ(result, test_ibaq_pep, "GFP_vs_WT", Rbbp4))
  expect_error(iBAQ(result[-(1)], test_ibaq_pep, "GFP_vs_WT", "Rbbp4"))
  expect_error(iBAQ(result, test_ibaq_pep[,-(6)], "GFP_vs_WT", "Rbbp4"))
  expect_error(iBAQ(result, test_ibaq_pep[,-(14)], "GFP_vs_WT", "Rbbp4"))
  expect_error(iBAQ(result, test_ibaq_pep, "test", "Rbbp4"))
  expect_error(iBAQ(result, test_ibaq_pep, "GFP_vs_WT", "test"))

  result_error <- result
  result_error$data <- result_error$data[,-(15:20)]
  expect_error(iBAQ(result_error, test_ibaq_pep, "GFP_vs_WT", "Rbbp4"))

  result_error2 <- result
  rowData(result_error2$sign) <- rowData(result_error2$sign)[,-(1)]
  expect_error(iBAQ(result_error2, test_ibaq_pep, "GFP_vs_WT", "Rbbp4"))

  result_error3 <- result
  rowData(result_error3$sign) <- rowData(result_error3$sign)[,-c(27:28,30)]
  expect_error(iBAQ(result_error3, test_ibaq_pep, "GFP_vs_WT", "Rbbp4"))
})

test_that("iBAQ returns a data.frame", {
  result <- LFQ(test_ibaq[,-(31:32)], GFPip_ExpDesign, "MinProb", "WT", "control", filter = c("Reverse", "Contaminant"),
                name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1)
  expect_is(iBAQ(result, test_ibaq_pep, "GFP_vs_WT", "Rbbp4", level = 1), "data.frame")
  expect_is(iBAQ(result, test_ibaq_pep, "GFP_vs_WT", "Rbbp4", level = 1L), "data.frame")
  expect_is(iBAQ(result, tibble::as_tibble(test_ibaq_pep), "GFP_vs_WT", "Rbbp4", level = 1), "data.frame")
})

