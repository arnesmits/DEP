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
  data %<>% unique_names(., "gene_name", "protein_id", delim = "[|]") %>% make_se(., cols, expdesign)
  filt <- filter_missval(data)
  norm <- norm_vsn(filt)
  imp <- imputation(norm, fun)
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
#' @param expdesign Data.frame, The experimental design object.
#' @param fun Character, Function used for data imputation based on \code{\link{impute}}.
#' @param control Character, The sample name to which the contrasts are generated (the control sample would be most appropriate).
#' @param type "all" or "control" The type of contrasts that will be generated.
#' @param filter Character, Name(s) of the column(s) to be filtered on.
#' @param name Character, Name of the column representing gene names.
#' @param ids "Character, Name of the column representing protein IDs.
#' @param alpha Numeric, sets the false discovery rate threshold.
#' @param lfc Numeric, sets the log fold change threshold.
#' @return A list of two data.frames: 1) \code{results} object containing the significant proteins, 2) \code{data} object containing the full dataset. Additionally, a \code{\link{rmarkdown}} report is generated and saved.
#' @examples
#' data <- UbIA_MS
#' expdesign <- ExpDesign_UbIA_MS
#' results <- LFQ(data, expdesign, "MinProb", "Ctrl", "control")
#' @export
LFQ <- function(data, expdesign, fun, control, type, filter = c("Reverse", "Potential.contaminant"), name = "Gene.names", ids = "Protein.IDs", alpha = 0.05, lfc = 1) {
  cols <- grep("^LFQ", colnames(data))
  cols_filt <- grep(paste("^", filter, "$", sep = "", collapse = "|"), colnames(data))

  if (!is.null(cols_filt)) {
    if (length(cols_filt) == 1) {
      data %<>% filter(.[,cols_filt] != "+")
    } else {
      data %<>% filter(!apply(.[,cols_filt] == "+", 1, any))
    }
  }

  data %<>% unique_names(., name, ids, delim = ";") %>% make_se(., cols, expdesign)
  filt <- filter_missval(data)
  norm <- norm_vsn(filt)
  imp <- imputation(norm, fun)
  lm <- linear_model(imp, control, type)
  sign <- cutoffs(lm, alpha, lfc)
  res <- results(sign)

  param <- data.frame(alpha, lfc)
  return(list(se = data, filt = filt, norm = norm, imputed = imp, lm = lm, sign = sign, results = res, param = param))
}

#' Generate a markdown report
#'
#' \code{report} is a wrapper function running the entire analysis workflow for label free quantification (LFQ)-based proteomics data.
#'
#' @param results List of SummerizedExperiments obtained from \code{\link{LFQ}} or \code{\link{TMT}} functions.
#' @return A \code{\link{rmarkdown}} report is generated and saved.
#' @examples
#' data <- UbIA_MS
#' expdesign <- ExpDesign_UbIA_MS
#'
#' results <- LFQ(data, expdesign, "MinProb", "Ctrl", "control")
#' report(results)
#'
#' @export
report <- function(results) {
  data <- results$se
  filt <- results$filt
  norm <- results$norm
  sign <- results$sign
  param <- results$param
  table <- results$results

  wd <- paste(getwd(), "/Report", sep = "")
  dir.create(wd)
  file <- paste(system.file(package = "proteomeR"), "/Report.Rmd", sep = "")
  rmarkdown::render(file, output_format = "all", output_dir = wd, quiet = T)

  write.table(table, paste(wd, "results.txt", sep = "/"), row.names = FALSE, sep = "\t")
}
