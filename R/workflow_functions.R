#' TMT workflow
#'
#' \code{TMT} is a wrapper function running the entire analysis workflow for TMT-based proteomics data.
#'
#' @param data Data.frame, The data object.
#' @param expdesign Data.frame, The experimental design object.
#' @param fun Character, Function used for data imputation based on \code{\link{impute}}.
#' @param control Character, The sample name to which the contrasts are generated (the control sample would be most appropriate).
#' @param type "all" or "control" The type of contrasts that will be generated.
#' @param alpha Numeric, sets the false discovery rate threshold.
#' @param lfc Numeric, sets the log fold change threshold.
#' @return A list of two data.frames: 1) \code{results} object containing the significant proteins, 2) \code{data} object containing the full dataset. Additionally, a \code{\link{rmarkdown}} report is generated and saved.
#' @examples
#'
#' TMT_res <- TMT()
#' @export
TMT <- function(data, expdesign, fun, control, type, alpha = 0.05, lfc = 1) {
  data <- data[-grep("###", data$gene_name),]
  cols <- grep("signal_sum", colnames(data))
  data %<>% unique_names(., "gene_name", "protein_id", delim = "[|]") %>% into_exprset_expdesign(., cols, expdesign)
  filt <- miss_val_filter(data)
  norm <- norm_vsn(filt)
  imp <- imputation_MSn(norm, fun)
  lm <- linear_model(imp, control, type)
  sign <- cutoffs(lm, alpha, lfc)
  res <- results(sign)

  wd <- paste(getwd(), "/Report", sep = "")
  dir.create(wd)
  file <- paste(system.file(package = "proteomeR"), "/Report.Rmd", sep = "")
  rmarkdown::render(file, output_format = "all", output_dir = wd, quiet = T)
  return(list(results = res, data = sign))
}

#' LFQ workflow
#'
#' \code{LFQ} is a wrapper function running the entire analysis workflow for label free quantification (LFQ)-based proteomics data.
#'
#' @param data Data.frame, The data object.
#' @param fun Character, Function used for data imputation based on \code{\link{impute}}.
#' @param control Character, The sample name to which the contrasts are generated (the control sample would be most appropriate).
#' @param type "all" or "control" The type of contrasts that will be generated.
#' @param alpha Numeric, sets the false discovery rate threshold.
#' @param lfc Numeric, sets the log fold change threshold.
#' @return A list of two data.frames: 1) \code{results} object containing the significant proteins, 2) \code{data} object containing the full dataset. Additionally, a \code{\link{rmarkdown}} report is generated and saved.
#' @examples
#' example <- UbIA_MS
#' example_results <- LFQ(example, "QRILC", "Con_", "control", alpha = 0.05, lfc = 1)
#' @export
LFQ <- function(data, fun, control, type, alpha = 0.05, lfc = 1) {
  data %<>% filter(Reverse != "+", Potential.contaminant != "+")
  cols <- grep("LFQ.", colnames(data))
  data %<>% unique_names(., "Gene.names", "Protein.IDs", delim = ";") %>% into_exprset(., cols)
  filt <- miss_val_filter(data)
  norm <- norm_vsn(filt)
  imp <- imputation_MSn(norm, fun)
  lm <- linear_model(imp, control, type)
  sign <- cutoffs(lm, alpha, lfc)
  res <- results(sign)

  wd <- paste(getwd(), "/Report", sep = "")
  dir.create(wd)
  file <- paste(system.file(package = "proteomeR"), "/Report.Rmd", sep = "")
  rmarkdown::render(file, output_format = "all", output_dir = wd, quiet = T)
  return(list(results = res, data = sign))
}
