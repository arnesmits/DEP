#' TMT workflow
#'
#' \code{TMT} is a wrapper function running the entire analysis workflow for TMT-based proteomics data.
#'
#' @param data Data.frame, The data object.
#' @param expdesign Data.frame, The experimental design object.
#' @param fun Character, Function used for data imputation based on \code{\link[MSnbase]{impute}}.
#' @param control Character, The sample name to which the contrasts are generated (the control sample would be most appropriate).
#' @param type "all" or "control" The type of contrasts that will be generated.
#' @param name Character, Name of the column representing gene names.
#' @param ids "Character, Name of the column representing protein IDs.
#' @param alpha Numeric, sets the false discovery rate threshold.
#' @param lfc Numeric, sets the log fold change threshold.
#' @return A list of two data.frames: 1) \code{results} object containing the significant proteins, 2) \code{data} object containing the full dataset.
#' @examples
#' \dontrun{
#'
#' TMT_res <- TMT()
#'
#' }
#' @export
TMT <- function(data, expdesign, fun, control, type, name = "gene_name", ids = "protein_id", alpha = 0.05, lfc = 1) {
  if(is.integer(alpha)) alpha <- as.numeric(alpha)
  if(is.integer(lfc)) lfc <- as.numeric(lfc)
  assert_that(is.data.frame(data), is.data.frame(expdesign), is.character(fun), is.character(control), is.character(type),
              is.character(name), is.character(ids), is.numeric(alpha), is.numeric(lfc))

  # Filter the data for Reverse hits (indicated by "###" in the gene_name)
  data <- data[-grep("###", data$gene_name),]
  # Make unique names and turn the data into a SummarizedExperiment
  cols <- grep("signal_sum", colnames(data)) # The reporter signal columns
  data %<>% unique_names(., name, ids, delim = "[|]") %>% make_se(., cols, expdesign)
  # Filter on missing values
  filt <- filter_missval(data)
  # Variance stabilization
  norm <- norm_vsn(filt)
  # Impute missing values
  imp <- imputation(norm, fun)
  # Test for differential expression by empirical Bayes moderation of a linear model and defined contrasts
  lm <- linear_model(imp, control, type)
  # Denote differential expressed proteins
  sign <- cutoffs(lm, alpha, lfc)
  # Generate a results table
  res <- results(sign)

  param <- data.frame(alpha, lfc)
  return(list(se = data, filt = filt, norm = norm, imputed = imp, lm = lm, sign = sign, results = res, param = param))
}

#' LFQ workflow
#'
#' \code{LFQ} is a wrapper function running the entire analysis workflow for label free quantification (LFQ)-based proteomics data.
#'
#' @param data Data.frame, The data object.
#' @param expdesign Data.frame, The experimental design object.
#' @param fun Character, Function used for data imputation based on \code{\link[MSnbase]{impute}}.
#' @param control Character, The sample name to which the contrasts are generated (the control sample would be most appropriate).
#' @param type "all" or "control" The type of contrasts that will be generated.
#' @param filter Character, Name(s) of the column(s) to be filtered on.
#' @param name Character, Name of the column representing gene names.
#' @param ids "Character, Name of the column representing protein IDs.
#' @param alpha Numeric, sets the false discovery rate threshold.
#' @param lfc Numeric, sets the log fold change threshold.
#' @return A list of two data.frames: 1) \code{results} object containing the significant proteins, 2) \code{data} object containing the full dataset.
#' @examples
#' \dontrun{
#'
#' data <- UbiLength
#' expdesign <- UbiLength_ExpDesign
#' results <- LFQ(data, expdesign, "MinProb", "Ctrl", "control")
#'
#' }
#' @export
LFQ <- function(data, expdesign, fun, control, type, filter = c("Reverse", "Potential.contaminant"), name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1) {
  if(is.integer(alpha)) alpha <- as.numeric(alpha)
  if(is.integer(lfc)) lfc <- as.numeric(lfc)
  assert_that(is.data.frame(data), is.data.frame(expdesign), is.character(fun), is.character(control), is.character(type),
              is.character(filter), is.character(name), is.character(ids), is.numeric(alpha), is.numeric(lfc))
  # Filter out the positive proteins (indicated by "+") in the pre-defined columns
  cols_filt <- grep(paste("^", filter, "$", sep = "", collapse = "|"), colnames(data)) # The columns to filter on
  if (!is.null(cols_filt)) {
    if (length(cols_filt) == 1) {
      data %<>% filter(.[,cols_filt] != "+")
    } else {
      data %<>% filter(!apply(.[,cols_filt] == "+", 1, any))
    }
  }

  # Make unique names and turn the data into a SummarizedExperiment
  cols <- grep("^LFQ", colnames(data)) # The LFQ intensity columns
  data %<>% unique_names(., name, ids, delim = ";") %>% make_se(., cols, expdesign)
  # Filter on missing values
  filt <- filter_missval(data)
  # Variance stabilization
  norm <- norm_vsn(filt)
  # Impute missing values
  imp <- imputation(norm, fun)
  # Test for differential expression by empirical Bayes moderation of a linear model and defined contrasts
  lm <- linear_model(imp, control, type)
  # Denote differential expressed proteins
  sign <- cutoffs(lm, alpha, lfc)
  # Generate a results table
  res <- results(sign)

  param <- data.frame(alpha, lfc)
  return(list(se = data, filt = filt, norm = norm, imputed = imp, lm = lm, sign = sign, results = res, param = param))
}

#' Generate a markdown report
#'
#' \code{report} is a wrapper function running the entire analysis workflow for label free quantification (LFQ)-based proteomics data.
#'
#' @param results List of SummerizedExperiments obtained from \code{\link{LFQ}} or \code{\link{TMT}} functions.
#' @return A \code{\link[rmarkdown]{rmarkdown}} report is generated and saved.
#' @examples
#' \dontrun{
#'
#' data <- UbiLength
#' expdesign <- UbiLength_ExpDesign
#'
#' results <- LFQ(data, expdesign, "MinProb", "Ctrl", "control")
#' report(results)
#'
#' }
#' @export
report <- function(results) {
  assert_that(is.list(results))

  # Extract the objects used in the Rmarkdown report from the results object
  data <- results$se
  filt <- results$filt
  norm <- results$norm
  sign <- results$sign
  param <- results$param
  table <- results$results

  # Render the Rmarkdown report
  wd <- paste(getwd(), "/Report", sep = "")
  dir.create(wd)
  file <- paste(system.file(package = "proteomeR"), "/Report.Rmd", sep = "")
  rmarkdown::render(file, output_format = "all", output_dir = wd, quiet = T)

  # Save the results table in a tab-delimited txt file
  write.table(table, paste(wd, "results.txt", sep = "/"), row.names = FALSE, sep = "\t")
}
