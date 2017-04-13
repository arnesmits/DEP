#' Unique names
#'
#' \code{unique_names} generates unique identifiers for a dataset based on a "name" and "id" column.
#'
#' @param data Data.frame, The data object for which unique names will be created.
#' @param name Character, Name of the column containing feature names.
#' @param ids Character, Name of the column containing feature IDs.
#' @param delim Character, Delimiter seperating the proteins within on protein group.
#' @return A data.frame containing an addiontal variable "name" that contains unique names.
#' @examples
#' example <- UbiLength
#' example_unique <- unique_names(example, "Gene.names", "Protein.IDs", delim = ";")
#' @export
unique_names <- function(data, name, ids, delim = ";") {
  # Show error if inputs are not the required classes
  assertthat::assert_that(is.data.frame(data), is.character(name), is.character(ids), is.character(delim))

  # Show error if inputs do not contain required columns
  if (length(grep(paste("^", name, "$", sep = ""), colnames(data))) < 1) {
    stop("name input is not present in data", call. = FALSE)
  }
  if (length(grep(paste("^", ids, "$", sep = ""), colnames(data))) < 1) {
    stop("ids input is not present in data", call. = FALSE)
  }

  # Select the name and id columns, take the first identifier per row and make unique names. If there is no name, the ID will be taken.
  names <- data %>% select(matches(paste("^", name, "$", sep = "")), matches(paste("^", ids, "$", sep = ""))) %>%
    mutate(name = gsub(paste(delim, ".*", sep = ""), "", .[,1]), ID = gsub(paste(delim, ".*", sep = ""), "", .[,2]), name = make.unique(ifelse(name == "", ID, name)))
  data <- left_join(data, names)
  return(data)
}

#' Data.frame to SummarizedExperiment parsing from column names
#'
#' \code{make_se_parse} creates a SummarizedExperiment object based on a single data.frame
#'
#' @param data Data.frame, Data object which will be turned into a SummarizedExperiment
#' @param columns Vector of integers, Column numbers that contain the assay data.
#' @return A SummarizedExperiment object.
#' @examples
#' example <- UbiLength %>% filter(Reverse != "+", Potential.contaminant != "+")
#' example_unique <- unique_names(example, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(example_unique))
#' example_se <- make_se_parse(example_unique, columns)
#' @export
make_se_parse <- function(data, columns) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(is.data.frame(data), is.integer(columns))

  # Show error if inputs do not contain required columns
  if (any(!c("name", "ID") %in% colnames(data))) {
    stop("'name' and/or 'ID' columns are not present in data;\nRun data_unique() to obtain the required data", call. = FALSE)
  }
  if (any(!apply(data[,columns], 2, is.numeric))) {
    stop("specified columns should be numeric;\nRun make_se_parse() with the appropriate columns as argument")
  }

  # Select the assay data
  rownames(data) <- data$name
  raw <- data[,columns]
  raw[raw == 0] <- NA
  raw <- log2(raw)
  colnames(raw) %<>% gsub(Rlibstree::getCommonPrefix(colnames(raw)), "", .) %>% make.names()

  # Select the rowData
  row_data <- data[,-columns]
  rownames(row_data) <- row_data$name

  # Generate the colData
  col_data <- colnames(raw) %>% data.frame(label = ., stringsAsFactors = F) %>%
    mutate(condition = substr(label,1,nchar(label)-1), replicate = substr(label, nchar(label), nchar(label))) %>%
    unite(ID, condition, replicate, remove = F)
  rownames(col_data) <- col_data$ID
  colnames(raw) <- col_data$ID[lapply(col_data$label, function(x) grep(x, colnames(raw))) %>% unlist()]
  raw <- raw[,!is.na(colnames(raw))]

  # Generate the SummarizedExperiment object
  dataset <- SummarizedExperiment(assays = as.matrix(raw), colData = col_data, rowData = row_data)
  return(dataset)
}

#' Data.frame to SummarizedExperiment using an experimental design
#'
#' \code{make_se} creates a SummarizedExperiment object based on two data.frames: the data and experimental design.
#'
#' @param data Data.frame, Data object which will be turned into a SummarizedExperiment
#' @param columns Vector of integers, Number of the columns that contain the assay data.
#' @param expdesign Data.frame, Experimental design object containing 'label', 'condition' and 'replicate' information
#' @return A SummarizedExperiment object.
#' @examples
#' example <- UbiLength %>% filter(Reverse != "+", Potential.contaminant != "+")
#' example_unique <- unique_names(example, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(example_unique))
#' exp_design <- UbiLength_ExpDesign
#' example_se <- make_se(example_unique, columns, exp_design)
#' @export
make_se <- function(data, columns, expdesign) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(is.data.frame(data), is.integer(columns), is.data.frame(expdesign))

  # Show error if inputs do not contain required columns
  if (any(!c("name", "ID") %in% colnames(data))) {
    stop("'name' and/or 'ID' columns are not present in data;\nRun data_unique() to obtain the required data", call. = FALSE)
  }
  if (any(!c("label", "condition", "replicate") %in% colnames(expdesign))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in expdesign", call. = FALSE)
  }
  if (any(!apply(data[,columns], 2, is.numeric))) {
    stop("specified columns should be numeric;\nRun make_se_parse() with the appropriate columns as argument")
  }

  # Select the assay data
  rownames(data) <- data$name
  raw <- data[,columns]
  raw[raw == 0] <- NA
  raw <- log2(raw)

  # Generate the colData from the experimental design and match these with the assay data
  expdesign %<>% mutate(condition = make.names(condition)) %>% unite(ID, condition, replicate, remove = F)
  rownames(expdesign) <- expdesign$ID
  colnames(raw) <- expdesign$ID[lapply(expdesign$label, function(x) grep(x, colnames(raw))) %>% unlist()]
  raw <- raw[,!is.na(colnames(raw))]

  # Select the rowData
  row_data <- data[,-columns]
  rownames(row_data) <- row_data$name

  # Generate the SummarizedExperiment object
  dataset <- SummarizedExperiment(assays = as.matrix(raw), colData = expdesign, rowData = row_data)
  return(dataset)
}

#' Filter on missing values
#'
#' \code{filter_missval} filters an SummarizedExperiment object based on missing values.
#'
#' @param data SummarizedExperiment, Data object which will be filtered.
#' @param thr Integer, sets the threshold for the allowed number of missing values per condition.
#' @return An filtered SummarizedExperiment object.
#' @examples
#' example <- UbiLength %>% filter(Reverse != "+", Potential.contaminant != "+")
#' example_unique <- unique_names(example, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(example_unique))
#' exp_design <- UbiLength_ExpDesign
#' example_se <- make_se(example_unique, columns, exp_design)
#'
#' example_stringent_filter <- filter_missval(example_se, thr = 0)
#' example_less_stringent_filter <- filter_missval(example_se, thr = 1)
#' @export
filter_missval <- function(data, thr = 0) {
  # Show error if inputs are not the required classes
  if(is.integer(thr)) thr <- as.numeric(thr)
  assertthat::assert_that(inherits(data, "SummarizedExperiment"), is.numeric(thr))

  # Show error if inputs do not contain required columns
  if (any(!c("name", "ID") %in% colnames(rowData(data)))) {
    stop("'name' and/or 'ID' columns are not present in data (rowData);\nRun data_unique() and make_se() (or make_se_parse()) to obtain the required data", call. = FALSE)
  }
  if (any(!c("label", "condition", "replicate") %in% colnames(colData(data)))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in data (colData);\nRun make_se() or make_se_parse() to obtain the required data ", call. = FALSE)
  }
  if (thr < 0 | thr > max(colData(data)$replicate)) {
    stop("invalid filter threshold applied;\nRun filter_missval() with a threshold ranging from 0 to the number of replicates")
  }

  # Make assay data binary (1 = valid value)
  bin_data <- assay(data)
  bin_data[!is.na(assay(data))] <- 1
  bin_data[is.na(assay(data))] <- 0

  # Filter data on the maximum allowed number of missing values per condition (defined by thr)
  keep <- bin_data %>% data.frame() %>% rownames_to_column(.) %>% gather(ID, value, 2:ncol(.)) %>%
    left_join(., data.frame(colData(data)), by = "ID") %>% group_by(rowname, condition) %>%
    summarize(miss_val = n()-sum(value)) %>% filter(miss_val <= thr) %>% spread(condition, miss_val)
  filt <- data[keep$rowname,]
  return(filt)
}

#' Data normalization using vsn
#'
#' \code{norm_vsn} normalizes an SummarizedExperiment object using \code{\link[vsn]{vsn-package}}.
#'
#' @param data SummarizedExperiment, Data object which will be normalized with log2-transformed assay data.
#' @return An normalized SummarizedExperiment object.
#' @examples
#' example <- UbiLength %>% filter(Reverse != "+", Potential.contaminant != "+")
#' example_unique <- unique_names(example, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(example_unique))
#' exp_design <- UbiLength_ExpDesign
#' example_se <- make_se(example_unique, columns, exp_design)
#'
#' example_filter <- filter_missval(example_se, thr = 0)
#' example_vsn <- norm_vsn(example_filter)
#' @export
norm_vsn <- function(data) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(data, "SummarizedExperiment"))

  # Variance stabilization transformation on assay data
  data_vsn <- data
  vsn.fit <- vsn::vsnMatrix(2^assay(data_vsn))
  assay(data_vsn) <- vsn::predict(vsn.fit, 2^assay(data_vsn))
  return(data_vsn)
}

#' Imputation by random draws from a manually defined distribution
#'
#' \code{manual_impute} missing values in a Summarized
#'
#' @param data SummarizedExperiment, Data object for which missing values will be imputed.
#' @param shift Integer
#' @param scale Integer
#' @return An imputed SummarizedExperiment object.
#' @examples
#' example <- UbiLength %>% filter(Reverse != "+", Potential.contaminant != "+")
#' example_unique <- unique_names(example, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(example_unique))
#' exp_design <- UbiLength_ExpDesign
#' example_se <- make_se(example_unique, columns, exp_design)
#'
#' example_filter <- filter_missval(example_se, thr = 0)
#' example_vsn <- norm_vsn(example_filter)
#'
#' example_impute_manual <- imputation(example_vsn, fun = "man", shift = 1.8, scale = 0.3)
#' @export
manual_impute <- function(data, scale = 0.3, shift = 1.8) {
  if (is.integer(scale)) scale <- is.numeric(scale)
  if (is.integer(shift)) shift <- is.numeric(shift)
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(data, "SummarizedExperiment"), is.numeric(scale), is.numeric(shift))

  # Get descriptive parameters of the current sample distributions
  stat <- assay(data) %>% data.frame() %>% rownames_to_column(.) %>% gather(samples, value, 2:ncol(.)) %>%
    filter(!is.na(value))  %>% group_by(samples) %>%
    summarise(mean = mean(value), median = median(value), sd = sd(value), n = n(), infin = nrow(assay(data))-n)
  # Impute missing values by random draws from a distribution
  # which is left-shifted by parameters 'shift' * sd and scaled by parameter 'scale' * sd.
  for(a in 1:nrow(stat)) {
    assay(data)[is.na(assay(data)[,stat$samples[a]]),stat$samples[a]] <-
      rnorm(stat$infin[a] , mean = stat$median[a] - shift * stat$sd[a], sd = stat$sd[a] * scale)
  }
  return(data)
}

#' Obtain a MSnSet object from a SummarizedExperiment object
#'
#' \code{se2msn} generats a MSnSet object from a SummarizedExperiment object
#'
#' @param data SummarizedExperiment, Data object which will be turned into a MSnSet object.
#' @return A MSnSet object.
#' @examples
#' example <- UbiLength %>% filter(Reverse != "+", Potential.contaminant != "+")
#' example_unique <- unique_names(example, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(example_unique))
#' exp_design <- UbiLength_ExpDesign
#' example_se <- make_se(example_unique, columns, exp_design)
#'
#' example_MSnSet <- se2msn(example_se)
#' @export
se2msn <- function(data) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(data, "SummarizedExperiment"))

  # Extract expression, feature and pheno data
  raw <- assay(data)
  feat_data <- data.frame(rowData(data))
  rownames(feat_data) <- feat_data$name
  pheno_data <- data.frame(colData(data))

  # Generate MSnSet object
  msn <- MSnbase::MSnSet(exprs = as.matrix(raw),
                         pData = Biobase::AnnotatedDataFrame(pheno_data),
                         fData = Biobase::AnnotatedDataFrame(feat_data))
  return(msn)
}

#' Impute missing values
#'
#' \code{imputation} imputes missing values based on \code{\link[MSnbase]{impute}}.
#'
#' @param data SummarizedExperiment, Data object for which missing values will be imputed.
#' @param fun "man", "QRILC", "MinDet", "MinProb", "min", "zero", "MLE", "bpca" or "knn", Function used for data imputation based on \code{\link[MSnbase]{impute}}.
#' @param ... Additional arguments for imputation functions as depicted in \code{\link{manual_impute}} \code{\link[MSnbase]{impute}}.
#' @return An imputed SummarizedExperiment object.
#' @examples
#' example <- UbiLength %>% filter(Reverse != "+", Potential.contaminant != "+")
#' example_unique <- unique_names(example, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(example_unique))
#' exp_design <- UbiLength_ExpDesign
#' example_se <- make_se(example_unique, columns, exp_design)
#'
#' example_filter <- filter_missval(example_se, thr = 0)
#' example_vsn <- norm_vsn(example_filter)
#'
#' example_impute_MinProb <- imputation(example_vsn, fun = "MinProb", q = 0.05)
#' example_impute_QRILC <- imputation(example_vsn, fun = "QRILC")
#'
#' example_impute_knn <- imputation(example_vsn, fun = "knn", k = 10, rowmax = 0.9)
#' example_impute_MLE <- imputation(example_vsn, fun = "MLE")
#'
#' example_impute_manual <- imputation(example_vsn, fun = "man", shift = 1.8, scale = 0.3)
#' @export
imputation <- function(data, fun, ...) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(data, "SummarizedExperiment"), is.character(fun))

  # Show error if inputs do not contain required columns
  if (any(!c("name", "ID") %in% colnames(rowData(data)))) {
    stop("'name' and/or 'ID' columns are not present in data (rowData);\nRun data_unique() and make_se() (or make_se_parse()) to obtain the required data", call. = FALSE)
  }
  if (!fun %in% c("man", MSnbase::imputeMethods())) {
    stop(paste("run imputation() with a valid function;\nValid functions are ", paste(c("man", MSnbase::imputeMethods()), collapse = "', '"), "", sep = "'"))
  }

  # if the "man" function is selected, use the manual impution method
  if (fun == "man") {
    imputed <- manual_impute(data, ...)
  }
  # else use the MSnSet::impute function
  else {
    MSnSet_data <- se2msn(data)
    MSnSet_imputed <- MSnbase::impute(MSnSet_data, method = fun, ...)
    assay(data) <- MSnbase::exprs(MSnSet_imputed)
    return(data)
  }
}

#' Linear model test
#'
#' \code{linear_model} performs a differential expression test based on linear models and empherical Bayes statistics (\code{\link[limma]{limma}}).
#'
#' @param data SummarizedExperiment, Data object for which the variance will be analyzed.
#' @param control Character, The condition to which the contrasts are generated (the control would be most appropriate).
#' @param type all" or "control" The type of contrasts that will be generated.
#' @return An SummarizedExperiment object containing FDR estimates of differential expression.
#' @examples
#' example <- UbiLength %>% filter(Reverse != "+", Potential.contaminant != "+")
#' example_unique <- unique_names(example, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(example_unique))
#' exp_design <- UbiLength_ExpDesign
#' example_se <- make_se(example_unique, columns, exp_design)
#'
#' example_filter <- filter_missval(example_se, thr = 0)
#' example_vsn <- norm_vsn(example_filter)
#' example_impute <- imputation(example_filter, fun = "MinProb", q = 0.01)
#'
#' example_lm <- linear_model(example_impute, "Ctrl", "control")
#' @export
linear_model <- function(data, control, type) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(data, "SummarizedExperiment"), is.character(control), is.character(type))

  # Show error if inputs do not contain required columns
  if (any(!c("name", "ID") %in% colnames(rowData(data)))) {
    stop("'name' and/or 'ID' columns are not present in data (rowData);\nRun data_unique() and make_se() (or make_se_parse()) to obtain the required data", call. = FALSE)
  }
  if (any(!c("label", "condition", "replicate") %in% colnames(colData(data)))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in data (colData);\nRun make_se() or make_se_parse() to obtain the required data ", call. = FALSE)
  }
  # Show error if inputs are not valid
  if (!type %in% c("all", "control")) {
    stop("Not a valid type, run linear_model() with a valid type\nValid types are: 'all', 'control'", call. = FALSE)
  }
  if (!control %in% unique(colData(data)$condition)) {
    stop("Not a valid control; Run linear_model() with a valid control", paste0("\nValid controls are: '", paste0(unique(colData(data)$condition), collapse = "', '"), "'"), call. = FALSE)
  }

  # Make an appropriate design matrix
  conditions <- factor(colData(data)$condition)
  design <- model.matrix(~ 0 + conditions)
  colnames(design) <- gsub("conditions", "", colnames(design))

  # Generate contrasts to be tested
  # Either make all possible combinations ("all") or only the contrast versus the control sample ("control")
  if(type == "all") {
    cntrst <- apply(combn(colnames(design), 2), 2, function(x) paste(x, collapse = " - "))
    # Make sure that contrast containing the control sample have the control as denominator
    flip <- grep(paste("^", control, sep = ""), cntrst)
    if(length(flip) >= 1) {
      cntrst[flip] %<>% gsub(paste(control, "- ", sep = " "), "", .) %>% paste(" - ", control, sep = "")
    }
  }
  if(type == "control") {
    cntrst <- paste(colnames(design)[!colnames(design) %in% control], control, sep = " - ")
  }
  # Print tested contrasts
  cat("Tested contrasts: \n")
  print(gsub(" - ", "_vs_", cntrst))

  # Test for differential expression by empirical Bayes moderation of a linear model on the predefined contrasts
  eB_fit <- eBayes(contrasts.fit(lmFit(assay(data), design = design), makeContrasts(contrasts = cntrst, levels = design)))

  # function to retrieve the results of the differential expression test using 'fdrtool'
  retrieve_fun <- function(comp, fit = eB_fit){
    res <- topTable(fit, sort.by = "t", coef = comp, number = Inf)
    fdr_res <- fdrtool::fdrtool(res$t, plot = FALSE, verbose = FALSE)
    res$qval <- fdr_res$qval
    res$lfdr <- fdr_res$lfdr
    res$comparison <- rep(comp, dim(res)[1])
    res <- rownames_to_column(res)
    return(res)
  }

  # Retrieve the differential expression test restuls
  limma_res <- map_df(cntrst, retrieve_fun)

  # Select the logFC and qval variables
  limma_res_small <- limma_res %>% select(rowname, logFC, qval, comparison) %>% mutate(comparison = gsub(" - ", "_vs_", comparison))
  # Obtain a wide table with all log2 fold changes
  table_diff <- limma_res_small %>% select(rowname, logFC, comparison) %>% spread(comparison, logFC)
  colnames(table_diff)[2:ncol(table_diff)] <- paste(colnames(table_diff)[2:ncol(table_diff)], "diff", sep = "_")
  # Obtain a wide table with all q-values
  table_padj <- limma_res_small %>% select(rowname, qval, comparison) %>% spread(comparison, qval)
  colnames(table_padj)[2:ncol(table_padj)] <- paste(colnames(table_padj)[2:ncol(table_padj)], "p.adj", sep = "_")
  # Join the two tables with the rowData
  table <- left_join(table_diff, table_padj, by = "rowname")
  rowData(data) <- merge(rowData(data), table, by.x = "name", by.y = "rowname")
  return(data)
}

#' Denote significant proteins
#'
#' \code{cutoffs} denotes significant proteins based on defined cutoffs.
#'
#' @param data SummarizedExperiment, Data object which will be filtered for significant proteins.
#' @param alpha Integer, sets the threshold for the false discovery rate (FDR).
#' @param lfc Integer, sets the threshold for the log fold change (lfc).
#' @return An SummarizedExperiment object containing a variable "sign" denoting the significant proteins.
#' @examples
#' example <- UbiLength %>% filter(Reverse != "+", Potential.contaminant != "+")
#' example_unique <- unique_names(example, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(example_unique))
#' exp_design <- UbiLength_ExpDesign
#' example_se <- make_se(example_unique, columns, exp_design)
#'
#' example_filter <- filter_missval(example_se, thr = 0)
#' example_vsn <- norm_vsn(example_filter)
#' example_impute <- imputation(example_filter, fun = "MinProb", q = 0.01)
#'
#' example_lm <- linear_model(example_impute, "Ctrl", "control")
#' example_sign <- cutoffs(example_lm, alpha = 0.05, lfc = 1)
#' significant_proteins <- example_sign[SummarizedExperiment::rowData(example_sign)$sign == "+", ]
#' nrow(significant_proteins)
#' @export
cutoffs <- function(data, alpha = 0.05, lfc = 1) {
  # Show error if inputs are not the required classes
  if(is.integer(alpha)) alpha <- as.numeric(alpha)
  if(is.integer(lfc)) lfc <- as.numeric(lfc)
  assertthat::assert_that(inherits(data, "SummarizedExperiment"), is.numeric(alpha), is.numeric(lfc))

  # Show error if inputs do not contain required columns
  if (any(!c("name", "ID") %in% colnames(rowData(data)))) {
    stop("'name' and/or 'ID' columns are not present in data", call. = FALSE)
  }
  if (length(grep("_p.adj|_diff", colnames(rowData(data)))) < 1) {
    stop("'[contrast]_diff' and/or '[contrast]_p.adj' columns are not present in data;\nRun linear_model() to obtain the required data", call. = FALSE)
  }

  row_data <- rowData(data)
  cols_p <- grep("_p.adj",colnames(row_data)) # get all columns with adjusted p-values
  cols_diff <- grep("_diff", colnames(row_data)) # get all columns with log2 fold changes

  # Denote differential expressed proteins by applying the alpha and logFC parameters per protein
  if(length(cols_p) == 1) {
    rowData(data)$sign <- ifelse(row_data[,cols_p] <= alpha & row_data[,cols_diff] >= lfc | row_data[,cols_p] <= alpha & row_data[,cols_diff] <= -lfc, "+", "")
    rowData(data)$contrast_sign <- rowData(data)$sign
    colnames(rowData(data))[ncol(rowData(data))] <- gsub("p.adj", "sign", colnames(row_data)[cols_p])
  }
  if(length(cols_p) > 1) {
    sign_df <- ifelse(row_data[,cols_p] %>% apply(., 2, function(x) x <= alpha) & row_data[,cols_diff] %>% apply(., 2, function(x) x >= lfc | x <= -lfc), "+", "")
    sign_df <- cbind(sign_df, sign = ifelse(apply(sign_df, 1, function(x) any(x == "+")),"+",""))
    colnames(sign_df) %<>% gsub("_p.adj","_sign",.)

    sign_df <- cbind(name = row_data$name, sign_df)
    rowData(data) <- merge(rowData(data), sign_df, by = "name")
  }
  return(data)
}

#' Generate a results table
#'
#' \code{results} generates a results table (data.frame) from a SummarizedExperiment object which has been generated by \code{\link{linear_model}} followed by \code{\link{cutoffs}}.
#'
#' @param data SummarizedExperiment, Data object which has been generated by \code{\link{linear_model}} and \code{\link{cutoffs}}.
#' @return An data.frame object containing all results variables from the performed analysis.
#' @examples
#' example <- UbiLength %>% filter(Reverse != "+", Potential.contaminant != "+")
#' example_unique <- unique_names(example, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(example_unique))
#' exp_design <- UbiLength_ExpDesign
#' example_se <- make_se(example_unique, columns, exp_design)
#'
#' example_filter <- filter_missval(example_se, thr = 0)
#' example_vsn <- norm_vsn(example_filter)
#' example_impute <- imputation(example_filter, fun = "MinProb", q = 0.01)
#'
#' example_lm <- linear_model(example_impute, "Ctrl", "control")
#' example_sign <- cutoffs(example_lm, alpha = 0.05, lfc = 1)
#' example_results <- results(example_sign)
#' glimpse(example_results)
#' @export
results <- function(data) {
  # Show error if inputs are not the required classes
  assert_that(inherits(data, "SummarizedExperiment"))

  # Show error if inputs do not contain required columns
  if (any(!c("name", "ID") %in% colnames(rowData(data)))) {
    stop("'name' and/or 'ID' columns are not present in data", call. = FALSE)
  }
  if (length(grep("_p.adj|_diff", colnames(rowData(data)))) < 1) {
    stop("'[contrast]_diff' and/or '[contrast]_p.adj' columns are not present in data;\nRun linear_model() to obtain the required data", call. = FALSE)
  }

  row_data <- data.frame(rowData(data))

  # Obtain average protein-centered enrichment values per condition
  rowData(data)$mean <- rowMeans(assay(data))
  centered <- assay(data) - rowData(data)$mean
  centered %<>%  data.frame(.) %>% rownames_to_column(.) %>% gather(ID, val, 2:ncol(.)) %>% left_join(., data.frame(colData(data)), by = "ID")
  centered %<>% group_by(rowname, condition) %>% summarize(val = mean(val)) %>% mutate(val = signif(val, digits = 3)) %>% spread(condition, val)
  colnames(centered)[2:ncol(centered)] %<>%  paste(., "_centered", sep = "")

  # Obtain average enrichments of conditions versus the control condition
  ratio <- row_data %>% column_to_rownames("name") %>% select(ends_with("diff")) %>% signif(., digits = 3) %>% rownames_to_column(.)
  colnames(ratio)[2:ncol(ratio)] %<>% gsub("_diff", "_ratio", .)
  df <- left_join(ratio, centered, by = "rowname")

  # Select the adjusted p-values and significance columns
  pval <- row_data %>% column_to_rownames("name") %>% select(ends_with("p.adj"), ends_with("sign")) %>% rownames_to_column(.)
  pval[,grep("p.adj", colnames(pval))] %<>%  signif(., digits = 3) %>% format(., scientific = T)

  # Join into a results table
  ids <- row_data %>% select(name, ID)
  table <- left_join(ids, pval, by = c("name" = "rowname"))
  table <- left_join(table, df, by = c("name" = "rowname"))
  return(table)
}
