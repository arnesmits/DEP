context("1 - Make unique names and generate a SummarizedExperiment")

test_that("make_unique throws error without valid input", {
  expect_error(make_unique("test_data", "Gene.names", "Protein.IDs"))
  expect_error(make_unique(test_data, "Gene.name", "Protein.IDs"))
  expect_error(make_unique(test_data, "Gene.names", "Protein.ID"))
})

test_that("make_unique returns a data.frame", {
  expect_is(make_unique(test_data, "Gene.names", "Protein.IDs"), "data.frame")
  expect_is(make_unique(tibble::as_tibble(test_data), "Gene.names", "Protein.IDs"), "data.frame")
})

test_that("make_unique returns unique names", {
  expect_false(any(duplicated(make_unique(test_data, "Gene.names", "Protein.IDs")$name)))
})

test_that("make_se_parse trows error without valid input", {
  expect_error(make_se_parse("test_unique", 21:32))
  expect_error(make_se_parse(test_unique, "21:32"))
  expect_error(make_se_parse(test_unique, 1:10))
  expect_error(make_se_parse(test_unique[,-(36:37)], 21:32))
})

test_that("make_se_parse returns a SummarizedExperiment", {
  expect_is(make_se_parse(test_unique, 21:32), "SummarizedExperiment")
  expect_is(make_se_parse(tibble::as_tibble(test_unique), 21:32), "SummarizedExperiment")
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
  expect_is(make_se(tibble::as_tibble(test_unique), 21:32, UbiLength_ExpDesign), "SummarizedExperiment")
  expect_is(make_se(test_unique, 21:32, tibble::as_tibble(UbiLength_ExpDesign)), "SummarizedExperiment")
})

