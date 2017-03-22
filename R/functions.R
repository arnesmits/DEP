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
#' example <- UbIA_MS
#' example_unique <- unique_names(example, "Gene.names", "Protein.IDs", delim = ";")
#' @export
unique_names <- function(data, name, ids, delim = ";") {
  names <- data %>% select(matches(name), matches(paste("^", ids, sep = ""))) %>%
    mutate(name = gsub(paste(delim, ".*", sep = ""), "", .[,1]), ID = gsub(paste(delim, ".*", sep = ""), "", .[,2]), name = make.unique(ifelse(name == "", ID, name)))
  data <- left_join(data, names)
  return(data)
}

#' Data.frame to SummerizedExperiment parsing from column names
#'
#' \code{make_se_parse} creates a SummerizedExperiment object based on a single data.frame
#'
#' @param data Data.frame, Data object which will be turned into a SummerizedExperiment
#' @param columns Vector of integers, Number of the columns that contain the assay data.
#' @return A SummerizedExperiment object.
#' @examples
#' example <- UbIA_MS %>% filter(Reverse != "+", Potential.contaminant != "+")
#' example_unique <- unique_names(example, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(example_unique))
#' example_exprset <- make_se_parse(example_unique, columns)
#' @export
make_se_parse <- function(data, columns) {
  rownames(data) <- data$name
  raw <- data[,columns]
  raw[raw == 0] <- NA
  raw <- log2(raw)

  colnames(raw) %<>% gsub(getCommonPrefix(colnames(raw)), "", .) %>% make.names()
  row_data <- data[,-columns]
  rownames(row_data) <- row_data$name

  col_data <- colnames(raw) %>% data.frame(ID = ., stringsAsFactors = F) %>% mutate(replicate = substr(ID, nchar(ID), nchar(ID)), condition = substr(ID,1,nchar(ID)-1))
  rownames(col_data) <- col_data$ID

  dataset <- SummarizedExperiment(assays = as.matrix(raw), colData = col_data, rowData = row_data)
  return(dataset)
}

#' Data.frame to SummerizedExperiment using an experimental design
#'
#' \code{make_se} creates a SummerizedExperiment object based on two data.frames: the data and experimental design.
#'
#' @param data Data.frame, Data object which will be turned into a SummerizedExperiment
#' @param columns Vector of integers, Number of the columns that contain the assay data.
#' @param expdesign Data.frame, Experimental design object
#' @return A SummerizedExperiment object.
#' @examples
#' example <- UbIA_MS %>% filter(Reverse != "+", Potential.contaminant != "+")
#' example_unique <- unique_names(example, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(example_unique))
#' exp_design <- ExpDesign_UbIA_MS
#' example_exprset <- make_se(example_unique, columns, exp_design)
#' @export
make_se <- function(data, columns, expdesign) {
  rownames(data) <- data$name
  raw <- data[,columns]
  raw[raw == 0] <- NA
  raw <- log2(raw)

  expdesign %<>% mutate(condition = make.names(condition)) %>% unite(ID, condition, replicate, remove = F)
  rownames(expdesign) <- expdesign$ID
  colnames(raw) <- expdesign$ID[lapply(expdesign$label, function(x) grep(x, colnames(raw))) %>% unlist()]
  raw <- raw[,!is.na(colnames(raw))]

  row_data <- data[,-columns]
  rownames(row_data) <- row_data$name
  dataset <- SummarizedExperiment(assays = as.matrix(raw), colData = expdesign, rowData = row_data)
  return(dataset)
}

#' Filter on missing values
#'
#' \code{filter_missval} filters an SummerizedExperiment object based on missing values.
#'
#' @param data SummerizedExperiment, Data object which will be filtered.
#' @param thr Integer, sets the threshold for the allowed number of missing values per condition.
#' @return An filtered SummerizedExperiment object.
#' @examples
#' example <- UbIA_MS %>% filter(Reverse != "+", Potential.contaminant != "+")
#' example_unique <- unique_names(example, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(example_unique))
#' exp_design <- ExpDesign_UbIA_MS
#' example_exprset <- make_se(example_unique, columns, exp_design)
#'
#' example_stringent_filter <- filter_missval(example_exprset, thr = 0)
#' example_less_stringent_filter <- filter_missval(example_exprset, thr = 1)
#' @export
filter_missval <- function(data, thr = 0) {
  bin_data <- assay(data)
  bin_data[!is.na(assay(data))] <- 1
  keep <- bin_data %>% data.frame() %>% rownames_to_column(.) %>% gather(ID, value, 2:ncol(.)) %>% left_join(., data.frame(colData(data)), by = "ID") %>% group_by(rowname, condition) %>%
    summarize(miss_val = n()-sum(value)) %>% filter(miss_val <= thr) %>% spread(condition, miss_val)
  filt <- data[keep$rowname,]
  return(filt)
}

#' Data normalization using vsn
#'
#' \code{norm_vsn} nomralizes an SummerizedExperiment object using \code{\link{vsn}}.
#'
#' @param data SummerizedExperiment, Data object which will be normalized.
#' @return An normalized SummerizedExperiment object.
#' @examples
#' example <- UbIA_MS %>% filter(Reverse != "+", Potential.contaminant != "+")
#' example_unique <- unique_names(example, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(example_unique))
#' exp_design <- ExpDesign_UbIA_MS
#' example_exprset <- make_se(example_unique, columns, exp_design)
#'
#' example_filter <- filter_missval(example_vsn, thr = 0)
#' example_vsn <- norm_vsn(example_filter)
#' @export
norm_vsn <- function(data) {
  data_vsn <- data
  vsn.fit <- vsnMatrix(2^assay(data_vsn))
  assay(data_vsn) <- predict(vsn.fit, 2^assay(data_vsn))
  return(data_vsn)
}

#' Impute missing values
#'
#' \code{imputation} imputes missing values based on a manual \code{\link{impute}}.
#'
#' @param data SummerizedExperiment, Data object for which missing values will be imputed.
#' @param fun "man", "QRILC", "MinDet", "MinProb", "min", "zero", "MLE", "bpca" or "knn", Function used for data imputation based on \code{\link{impute}}.
#' @param ... Additional arguments for "man" (scale and shift) or for other functions as depicted in \code{\link{impute}}.
#' @return An imputed SummerizedExperiment object.
#' @examples
#' example <- UbIA_MS %>% filter(Reverse != "+", Potential.contaminant != "+")
#' example_unique <- unique_names(example, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(example_unique))
#' exp_design <- ExpDesign_UbIA_MS
#' example_exprset <- make_se(example_unique, columns, exp_design)
#'
#' example_filter <- filter_missval(example_exprset, thr = 0)
#' example_vsn <- norm_vsn(example_filter)
#'
#' example_impute_MinProb <- imputation(example_filter, fun = "MinProb", q = 0.05)
#' example_impute_QRILC <- imputation(example_filter, fun = "QRILC")
#'
#' example_impute_knn <- imputation(example_filter, fun = "knn", k = 10, rowmax = 0.9)
#' example_impute_MLE <- imputation(example_filter, fun = "MLE")
#'
#' example_impute_manual <- imputation(example_filter, fun = "man", shift = 1.8, scale = 0.3)
#' @export
imputation <- function(data, fun, ...) {

  manual_impute <- function(data, scale = 0.3, shift = 1.8) {
    stat <- assay(data) %>% data.frame() %>% rownames_to_column(.) %>% gather(samples, value, 2:ncol(.)) %>% filter(!is.na(value))  %>% group_by(samples) %>%
      summarise(mean = mean(value), median = median(value), sd = sd(value), n = n(), infin = nrow(assay(data))-n)
    for(a in 1:nrow(stat)) {
      assay(data)[is.na(assay(data)[,stat$samples[a]]),stat$samples[a]] <- rnorm(stat$infin[a] , mean = stat$median[a] - shift * stat$sd[a], sd = stat$sd[a] * scale)
    }
    rowData(data)$mean <- rowMeans(assay(data))
    return(data)
  }
  se2msn <- function(data) {
    raw <- assay(data)
    feat_data <- data.frame(rowData(data))
    rownames(feat_data) <- feat_data$name
    pheno_data <- data.frame(colData(data))
    msn <- MSnSet(exprs = as.matrix(raw), pData = AnnotatedDataFrame(pheno_data), fData = AnnotatedDataFrame(feat_data))
    return(msn)
  }

  if (fun == "man") {
    imputed <- manual_impute(data, ...)
  }
  else {
    MSnSet_data <- se2msn(data)
    MSnSet_imputed <- impute(MSnSet_data, method = fun, ...)
    assay(data) <- exprs(MSnSet_imputed)
    rowData(data)$mean <- rowMeans(assay(data))
    return(data)
  }
}

#' Linear model test
#'
#' \code{linear_model} performs a differential expression test based on linear models and empherical Bayes statistics (\code{\link{limma}}).
#'
#' @param data SummerizedExperiment, Data object for which the variance will be analyzed.
#' @param control Character, The condition to which the contrasts are generated (the control would be most appropriate).
#' @param type all" or "control" The type of contrasts that will be generated.
#' @return An SummerizedExperiment object containing FDR estimates of differential expression.
#' @examples
#' example <- UbIA_MS %>% filter(Reverse != "+", Potential.contaminant != "+")
#' example_unique <- unique_names(example, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(example_unique))
#' exp_design <- ExpDesign_UbIA_MS
#' example_exprset <- make_se(example_unique, columns, exp_design)
#'
#' example_filter <- filter_missval(example_exprset, thr = 0)
#' example_vsn <- norm_vsn(example_filter)
#' example_impute <- imputation(example_filter, fun = "MinProb", q = 0.01)
#'
#' example_lm <- linear_model(example_impute, "Ctrl", "control")
#' @export
linear_model <- function(data, control, type) {
  conditions <- factor(colData(data)$condition)
  design <- model.matrix(~ 0 + conditions)
  colnames(design) <- gsub("conditions", "", colnames(design))

  if(type == "all") {
    cntrst <- apply(combn(colnames(design), 2), 2, function(x) paste(x, collapse = " - "))
    flip <- grep(paste("^", control, sep = ""), cntrst)
    if(length(flip) >= 1) {
      cntrst[flip] %<>% gsub(paste(control, "- ", sep = " "), "", .) %>% paste(" - ", control, sep = "")
    }
  }
  if(type == "control") {
    cntrst <- paste(colnames(design)[!colnames(design) %in% control], control, sep = " - ")
  }
  cat("Tested contrasts: \n")
  print(gsub(" - ", "_vs_", cntrst))

  eB_fit <- eBayes(contrasts.fit(lmFit(assay(data), design = design), makeContrasts(contrasts = cntrst, levels = design)))

  retrieve_fun <- function(comp, fit = eB_fit){
    res <- topTable(fit, sort.by = "t", coef = comp, number = Inf)
    fdr_res <- fdrtool(res$t, plot = FALSE, verbose = FALSE)
    res$qval <- fdr_res$qval
    res$lfdr <- fdr_res$lfdr
    res$comparison <- rep(comp, dim(res)[1])
    res <- rownames_to_column(res)
    return(res)
  }

  limma_res <- map_df(cntrst, retrieve_fun)

  limma_res_small <- limma_res %>% select(rowname, logFC, qval, comparison) %>% mutate(comparison = gsub(" - ", "_vs_", comparison))
  table_diff <- limma_res_small %>% select(rowname, logFC, comparison) %>% spread(comparison, logFC)
  colnames(table_diff)[2:ncol(table_diff)] <- paste(colnames(table_diff)[2:ncol(table_diff)], "diff", sep = "_")
  table_padj <- limma_res_small %>% select(rowname, qval, comparison) %>% spread(comparison, qval)
  colnames(table_padj)[2:ncol(table_padj)] <- paste(colnames(table_padj)[2:ncol(table_padj)], "p.adj", sep = "_")
  table <- left_join(table_diff, table_padj, by = "rowname")
  rowData(data) <- merge(rowData(data), table, by.x = "name", by.y = "rowname")
  return(data)
}

#' ANOVA and post-hoc Tukey
#'
#' \code{anova_tukey} performs an Analysis of Variance fit (\code{\link{aov}}) and subsequent Tukey Honest Significant Differences analysis (\code{\link{TukeyHSD}}).
#'
#' @param data SummerizedExperiment, Data object for which the variance will be analyzed.
#' @param control Character, The condition to which the contrasts are generated (the control would be most appropriate).
#' @param type all" or "control" The type of contrasts that will be generated.
#' @return An SummerizedExperiment object containing FDR estimates of differential expression.
#' @examples
#' example <- UbIA_MS %>% filter(Reverse != "+", Potential.contaminant != "+")
#' example_unique <- unique_names(example, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(example_unique))
#' exp_design <- ExpDesign_UbIA_MS
#' example_exprset <- make_se(example_unique, columns, exp_design)
#'
#' example_filter <- filter_missval(example_exprset, thr = 0)
#' example_vsn <- norm_vsn(example_filter)
#' example_impute <- imputation(example_filter, fun = "MinProb", q = 0.01)
#'
#' example_anova_tukey <- anova_tukey(example_impute, "Ctrl", "control")
#' @export
anova_tukey <- function(data, control, type) {
  long <- assay(data) %>% data.frame() %>% mutate(name = rownames(.)) %>% gather(ID, val, 1:(ncol(.)-1)) %>% left_join(., data.frame(colData(data)), by = "ID")

  cat("  ANOVA test \n")
  anova_p <- long %>% group_by(name) %>% do(anova = aov(val ~ condition, data = .) %>% summary(.) %>% .[[1]] %>% .$`Pr(>F)` %>% .[1])
  anova_p %<>% ungroup(name) %>% mutate(p = unlist(anova), padj = p.adjust(p, method = "BH")) %>% select(-anova) %>% data.frame()
  rowData(data) <- merge(rowData(data), anova_p, by = "name")

  cat("\n  Post-hoc test (Tukey) \n")
  tukey <- long %>% group_by(name) %>% do(tukey = aov(val ~ condition, data = .) %>% TukeyHSD(.) %>% .$condition %>% .[,c(1,4)])
  tukey_df <- unlist(tukey$tukey) %>% matrix(., nrow = nrow(data), byrow = T) %>% data.frame()
  colnames(tukey_df) <- c(paste(tukey$tukey[[1]] %>% row.names(),"diff",sep="_"),paste(tukey$tukey[[1]] %>% row.names(),"p.adj",sep="_"))
  tukey_df <- cbind(name = tukey$name, tukey_df)

  if(type == "control") {
    cols_diff <- grep(paste(control, ".*_diff", sep = ""), colnames(tukey_df))
    cols_p <- grep(paste(control, ".*_p.adj", sep = ""), colnames(tukey_df))
    tukey_df <- tukey_df[,c(1,cols_diff,cols_p)]

    flip_col <- grep(paste("^", control, ".*_diff", sep = ""), colnames(tukey_df))
    if(length(flip_col) >= 1) {
      tukey_df[,flip_col] <- -tukey_df[,flip_col]
      colnames(tukey_df)[flip_col] <- gsub(paste(control, "-", sep = ""), "", colnames(tukey_df)[flip_col]) %>% gsub("_diff", "", .) %>% paste("-", control, "_diff", sep = "")
    }
  }

  cols_p <- grep("p.adj",colnames(tukey_df))
  order <- tukey_df %>% select(cols_p) %>% names() %>% c("name",.)
  tukey_padj <- tukey_df[,c(1,cols_p)] %>% gather(., comparison, p, 2:ncol(.)) %>% mutate(padj = p.adjust(p, method = "BH")) %>% select(-p) %>% spread(., comparison, padj)
  final <- merge(tukey_df[,-cols_p],tukey_padj[,order], by = "name")
  rowData(data) <- merge(rowData(data), final, by="name")
  return(data)
}

#' Denote significant proteins
#'
#' \code{cutoffs} denotes significant proteins based on defined cutoffs.
#'
#' @param data SummerizedExperiment, Data object which will be filtered for significant proteins.
#' @param alpha Integer, sets the threshold for the false discovery rate (FDR).
#' @param lfc Integer, sets the threshold for the log fold change (lfc).
#' @return An SummerizedExperiment object containing a variable "sign" denoting the significant proteins.
#' @examples
#' example <- UbIA_MS %>% filter(Reverse != "+", Potential.contaminant != "+")
#' example_unique <- unique_names(example, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(example_unique))
#' exp_design <- ExpDesign_UbIA_MS
#' example_exprset <- make_se(example_unique, columns, exp_design)
#'
#' example_filter <- filter_missval(example_exprset, thr = 0)
#' example_vsn <- norm_vsn(example_filter)
#' example_impute <- imputation(example_filter, fun = "MinProb", q = 0.01)
#'
#' example_lm <- linear_model(example_impute, "Ctrl", "control")
#' example_sign <- cutoffs(example_lm, alpha = 0.05, lfc = 1)
#' significant_proteins <- example_sign[rowData(example_sign)$sign == "+", ]
#' nrow(significant_proteins)
#' @export
cutoffs <- function(data, alpha = 0.05, lfc = 1) {
  row_data <- rowData(data)
  cols_p <- grep("_p.adj",colnames(row_data))
  cols_diff <- grep("_diff", colnames(row_data))
  sign_df <- ifelse(row_data[,cols_p] %>% apply(., 2, function(x) x <= alpha) & row_data[,cols_diff] %>% apply(., 2, function(x) x >= lfc | x <= -lfc), "+", "")
  sign_df <- cbind(sign_df, sign = ifelse(apply(sign_df, 1, function(x) any(x == "+")),"+",""))
  colnames(sign_df) %<>% gsub("_p.adj","_sign",.)
  sign_df <- cbind(name = row_data$name, sign_df)
  rowData(data) <- merge(rowData(data), sign_df, by = "name")
  return(data)
}

#' Generate a results table
#'
#' \code{results} generates a results table (data.frame) from a SummerizedExperiment object which has been generated by \code{\link{linear_model}} followed by \code{\link{cutoffs}}.
#'
#' @param data SummerizedExperiment, Data object which has been generated by \code{\link{linear_model}} and \code{\link{cutoffs}}.
#' @return An data.frame object containing all results variables from the performed analysis.
#' @examples
#' example <- UbIA_MS %>% filter(Reverse != "+", Potential.contaminant != "+")
#' example_unique <- unique_names(example, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(example_unique))
#' exp_design <- ExpDesign_UbIA_MS
#' example_exprset <- make_se(example_unique, columns, exp_design)
#'
#' example_filter <- filter_missval(example_exprset, thr = 0)
#' example_vsn <- norm_vsn(example_filter)
#' example_impute <- imputation(example_filter, fun = "MinProb", q = 0.01)
#'
#' example_lm <- linear_model(example_impute, "Ctrl", "control")
#' example_sign <- cutoffs(example_lm, alpha = 0.05, lfc = 1)
#' example_results <- results(example_sign)
#' glimpse(example_results)
#' @export
results <- function(data) {
  row_data <- data.frame(rowData(data))
  centered <- assay(data) - rowData(data)$mean
  centered %<>%  data.frame(.) %>% rownames_to_column(.) %>% gather(ID, val, 2:ncol(.)) %>% left_join(., data.frame(colData(data)), by = "ID")
  centered %<>% group_by(rowname, condition) %>% summarize(val = mean(val)) %>% mutate(val = signif(val, digits = 3)) %>% spread(condition, val)
  colnames(centered)[2:ncol(centered)] %<>%  paste(., "_centered", sep = "")

  ratio <- row_data %>% column_to_rownames("name") %>% .[,grep("_diff$", colnames(.))] %>% signif(., digits = 3) %>% rownames_to_column(.)
  colnames(ratio)[2:ncol(ratio)] %<>% gsub("_diff", "_ratio", .)
  df <- left_join(ratio, centered, by = "rowname")

  ids <- row_data %>% select(name, ID)

  pval <- row_data %>% column_to_rownames("name") %>% .[,grep("p.adj$|sign$", colnames(.))] %>% rownames_to_column(.)
  pval[,grep("p.adj", colnames(pval))] %<>%  signif(., digits = 3) %>% format(., scientific = T)

  table <- left_join(ids, pval, by = c("name" = "rowname"))
  table <- left_join(table, df, by = c("name" = "rowname"))
  return(table)
}
