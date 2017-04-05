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
  # Show error if inputs are not the required classes
  if(is.integer(alpha)) alpha <- as.numeric(alpha)
  if(is.integer(lfc)) lfc <- as.numeric(lfc)
  assertthat::assert_that(is.data.frame(data), is.data.frame(expdesign), is.character(fun), is.character(control), is.character(type),
              is.character(name), is.character(ids), is.numeric(alpha), is.numeric(lfc))

  # Show error if inputs do not contain required columns
  if (length(grep(paste("^", name, "$", sep = ""), colnames(data))) < 1) {
    stop("name input is not present in data", call. = FALSE)
  }
  if (length(grep(paste("^", ids, "$", sep = ""), colnames(data))) < 1) {
    stop("ids input is not present in data", call. = FALSE)
  }
  if (any(!c("label", "condition", "replicate") %in% colnames(expdesign))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in expdesign", call. = FALSE)
  }
  # Show error if inputs are not valid
  if (!type %in% c("all", "control")) {
    stop("Not a valid type, run TMT() with a valid type\nValid types are: 'all', 'control'", call. = FALSE)
  }

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
  # Show error if inputs are not the required classes
  if(is.integer(alpha)) alpha <- as.numeric(alpha)
  if(is.integer(lfc)) lfc <- as.numeric(lfc)
  assertthat::assert_that(is.data.frame(data), is.data.frame(expdesign), is.character(fun), is.character(control), is.character(type),
              is.character(filter), is.character(name), is.character(ids), is.numeric(alpha), is.numeric(lfc))

  # Show error if inputs do not contain required columns
  if (length(grep(paste("^", name, "$", sep = ""), colnames(data))) < 1) {
    stop("name input is not present in data", call. = FALSE)
  }
  if (length(grep(paste("^", ids, "$", sep = ""), colnames(data))) < 1) {
    stop("ids input is not present in data", call. = FALSE)
  }
  if (length(grep(paste("^", filter, "$", sep = "", collapse = "|"), colnames(data))) < 1) {
    stop("filter input is not present in data", call. = FALSE)
  }
  if (any(!c("label", "condition", "replicate") %in% colnames(expdesign))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in expdesign", call. = FALSE)
  }
  # Show error if inputs are not valid
  if (!type %in% c("all", "control")) {
    stop("Not a valid type, run LFQ() with a valid type\nValid types are: 'all', 'control'", call. = FALSE)
  }

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
  data_unique <- unique_names(data, name, ids, delim = ";")
  se <- make_se(data_unique, cols, expdesign)
  # Filter on missing values
  filt <- filter_missval(se)
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
  return(list(data = data_unique, se = se, filt = filt, norm = norm, imputed = imp, lm = lm, sign = sign, results = res, param = param))
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
  # Show error if inputs are not the required classes
  assertthat::assert_that(is.list(results))

  # Show error in case that the required objects are not present in the list object 'results'
  if(any(!c("data", "se", "filt", "norm", "imputed", "lm", "sign", "results", "param") %in% names(results))) {
    stop("No valid input; Run report() with appropriate input generated by LFQ() or TMT()", call. = FALSE)
  }

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

#' iBAQ workflow
#'
#' \code{iBAQ} is a wrapper function running the entire analysis workflow for stoichiometry analysis using intensity-based absolute quantification (iBAQ)-based proteomics data.
#'
#' @param results List of SummerizedExperiments obtained from \code{\link{LFQ}} function.
#' @param peptides Data.frame, Peptide table from MaxQuant ("peptides.txt").
#' @param contrast Character, The specific contrast to calculate the stroichiometry for.
#' @param bait Character, The name of the protein to which all other proteins will be scaled for the relative stoichiometry
#' @param level Numerical, The level to which the bait will be scaled
#' @return A data.frame with the relative stoichiometry data
#' @examples
#' \dontrun{
#' # load data and test for differentially enriched proteins
#' data <- GFPip
#' expdesign <- GFPip_ExpDesign
#' results <- LFQ(data, expdesign, "MinProb", "WT", "control", filter = c("Reverse", "Contaminant"), alpha = 0.05, lfc = 4.5)
#'
#' # load peptide data and perform iBAQ-based stoichiometry analysis
#' peptides <- GFPip_pep
#' iBAQ(results, peptides, contrast = "GFP_vs_WT", bait = "Suz12")
#'
#' }
#' @export
iBAQ <- function(results, peptides, contrast, bait, level = 1) {
  # Show error if inputs are not the required classes
  if(is.integer(level)) level <- as.numeric(level)
  assertthat::assert_that(is.list(results), is.data.frame(peptides), is.character(contrast), is.character(bait), is.numeric(level))

  # Show error in case that the required objects are not present in the list object 'results'
  if(any(!c("data", "se", "filt", "norm", "imputed", "lm", "sign", "results", "param") %in% names(results))) {
    stop("No valid input; Run iBAQ() with appropriate input generated by LFQ()", call. = FALSE)
  }
  if (length(grep("iBAQ.", colnames(results$data))) < 1) {
    stop("'iBAQ' columns are not present in results", call. = FALSE)
  }
  if (any(!c("Protein.group.IDs", "Unique..Groups.") %in% colnames(peptides))) {
    stop("'Protein.group.IDs' and/or 'Unique..Groups.' columns are not present in peptides", call. = FALSE)
  }
  if (any(!c("name", "ID") %in% colnames(rowData(results$sign)))) {
    stop("'name' and/or 'ID' columns are not present in results$sign", call. = FALSE)
  }
  if (length(grep("_sign|_diff", colnames(rowData(results$sign)))) < 1) {
    stop("'[contrast]_sign' and/or '[contrast]_diff' columns are not present in results$sign", call. = FALSE)
  }

  # Merge iBAQ values for peptides with shared peptides
  ibaq <- merge_ibaq(results$data, peptides)

  # Calculate relative stoichiometry compared to the bait
  stoi <- stoichiometry(results$sign, ibaq, contrast, bait, level)

  # Plot stoichiometry
  p1 <- plot_stoi(stoi)
  print(p1)

  return(stoi)
}
