context("1 - Make unique names and generate a SummarizedExperiment")

test_that("unique_names throws error without valid input", {
  expect_error(unique_names("test_data", "Gene.names", "Protein.IDs"))
  expect_error(unique_names(test_data, "Gene.name", "Protein.IDs"))
  expect_error(unique_names(test_data, "Gene.names", "Protein.ID"))
})

test_that("unique_names returns a data.frame", {
  expect_is(unique_names(test_data, "Gene.names", "Protein.IDs"), "data.frame")
})

test_that("unique_names returns unique names", {
  expect_false(any(duplicated(unique_names(test_data, "Gene.names", "Protein.IDs")$name)))
})

test_that("make_se_parse trows error without valid input", {
  expect_error(make_se_parse("test_unique", 21:32))
  expect_error(make_se_parse(test_unique, "21:32"))
  expect_error(make_se_parse(test_unique, 1:10))
  expect_error(make_se_parse(test_unique[,-(36:37)], 21:32))
})

test_that("make_se_parse returns a SummarizedExperiment", {
  expect_is(make_se_parse(test_unique, 21:32), "SummarizedExperiment")
})

test_that("make_se trows error without valid input", {
  expect_error(make_se("test_unique", 21:32, UbiLength_ExpDesign))
  expect_error(make_se(test_unique, "21:32", UbiLength_ExpDesign))
  expect_error(make_se("test_unique", 21:32, "UbiLength_ExpDesign"))

  expect_error(make_se(test_unique, 1:10, UbiLength_ExpDesign))
  expect_error(make_se(test_unique[,-(36:37)], 21:32, UbiLength_ExpDesign))
  expect_error(make_se(test_unique, 21:32, UbiLength_ExpDesign[,-(2)]))
})

test_that("make_se returns a SummarizedExperiment", {
  expect_is(make_se(test_unique, 21:32, UbiLength_ExpDesign), "SummarizedExperiment")
})

