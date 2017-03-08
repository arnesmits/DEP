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

#' Data.frame to ExpressionSet
#'
#' \code{into_exprset} creates an ExpressionSet object based on a data.frame.
#'
#' @param data Data.frame, The data object which will be turned into an ExpressionSet.
#' @param columns Vector of integers, Number of the columns that contain the Expression data.
#' @return An ExpressionSet object.
#' @examples
#' example <- UbIA_MS %>% filter(Reverse != "+", Potential.contaminant != "+")
#' example_unique <- unique_names(example, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(example_unique))
#' example_exprset <- into_exprset(example_unique, columns)
#' @export
into_exprset <- function(data, columns) {
  rownames(data) <- data$name
  raw <- data[,columns]
  raw[raw == 0] <- NA
  raw <- log2(raw)
  colnames(raw) %<>% gsub(getCommonPrefix(colnames(raw)), "", .) %>% make.names()
  feature_data <- data[,-columns]
  rownames(feature_data) <- feature_data$name
  pheno_data <- colnames(raw) %>% data.frame(ID = ., stringsAsFactors = F) %>% mutate(replicate = substr(ID, nchar(ID), nchar(ID)), sample = substr(ID,1,nchar(ID)-1))
  rownames(pheno_data) <- pheno_data$ID
  dataset <- ExpressionSet(as.matrix(raw), phenoData = AnnotatedDataFrame(pheno_data), featureData = AnnotatedDataFrame(feature_data))
  return(dataset)
}

#' Data.frame to ExpressionSet using an experimental design
#'
#' \code{into_exprset_expdesign} creates an ExpressionSet object based on two data.frames: the data and experimental design.
#'
#' @param data Data.frame, Data object which will be turned into an ExpressionSet.
#' @param columns Vector of integers, Number of the columns that contain the Expression data.
#' @param expdesign Data.frame, Experimental design object
#' @return An ExpressionSet object.
#' @examples
#' example <- UbIA_MS %>% filter(Reverse != "+", Potential.contaminant != "+")
#' example_unique <- unique_names(example, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(example_unique))
#' exp_design <- ExpDesign_UbIA_MS
#' example_exprset <- into_exprset_expdesign(example_unique, columns, exp_design)
#' @export
into_exprset_expdesign <- function(data, columns, expdesign) {
  rownames(data) <- data$name
  raw <- data[,columns]
  raw[raw == 0] <- NA
  raw <- log2(raw)

  expdesign %<>% mutate(sample = make.names(sample)) %>% unite(ID, sample, replicate, remove = F)
  rownames(expdesign) <- expdesign$ID
  colnames(raw) <- expdesign$ID[lapply(expdesign$label, function(x) grep(x, colnames(raw))) %>% unlist()]
  raw <- raw[,!is.na(colnames(raw))]

  feature_data <- data[,-columns]
  rownames(feature_data) <- feature_data$name
  dataset <- ExpressionSet(as.matrix(raw), phenoData = AnnotatedDataFrame(expdesign), featureData = AnnotatedDataFrame(feature_data))
  return(dataset)
}

#' Data normalization using vsn
#'
#' \code{norm_vsn} nomralizes an ExpressionSet object using \code{\link{vsn}}.
#'
#' @param data ExpressionSet, Data object which will be normalized.
#' @return An normalized ExpressionSet object.
#' @examples
#' example <- UbIA_MS %>% filter(Reverse != "+", Potential.contaminant != "+")
#' example_unique <- unique_names(example, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(example_unique))
#' example_exprset <- into_exprset(example_unique, columns)
#'
#' example_vsn <- norm_vsn(example_exprset)
#' @export
norm_vsn <- function(data) {
  data_vsn <- data
  vsn.fit <- vsnMatrix(2^exprs(data_vsn))
  exprs(data_vsn) <- predict(vsn.fit, 2^exprs(data_vsn))
  return(data_vsn)
}

#' Filter on missing values
#'
#' \code{miss_val_filter} filters an ExpressionSet object based on missing values.
#'
#' @param data ExpressionSet, Data object which will be filtered.
#' @param thr Integer, sets the threshold for the allowed number of missing values per condition.
#' @return An filtered ExpressionSet object.
#' @examples
#' example <- UbIA_MS %>% filter(Reverse != "+", Potential.contaminant != "+")
#' example_unique <- unique_names(example, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(example_unique))
#' example_exprset <- into_exprset(example_unique, columns)
#'
#' example_vsn <- norm_vsn(example_exprset)
#' example_stringent_filter <- miss_val_filter(example_vsn, thr = 0)
#' example_less_stringent_filter <- miss_val_filter(example_vsn, thr = 1)
#' @export
miss_val_filter <- function(data, thr = 0) {
  bin_data <- exprs(data)
  bin_data[!is.na(exprs(data))] <- 1
  keep <- bin_data %>% data.frame() %>% mutate(name = rownames(.)) %>% gather(ID, value, 1:(ncol(.)-1)) %>% left_join(., pData(data), by = "ID") %>% group_by(name, sample) %>%
    summarize(miss_val = n()-sum(value)) %>% filter(miss_val <= thr) %>% spread(sample, miss_val)
  filt <- data[keep$name,]
  return(filt)
}

#' Impute missing values manually
#'
#' \code{imputation_perseus} imputes missing values by a distribution of values which is left-shifted and scaled from the original distribution.
#'
#' @param data ExpressionSet, Data object for which missing values will be imputed.
#' @param shift Integer, sets the left-shift based on the original distribution (1 corresponds to 1 sd left-shift from median).
#' @param scale Integer, sets the widht based on the original distribution sd.
#' @return An imputed ExpressionSet object.
#' @examples
#' example <- UbIA_MS %>% filter(Reverse != "+", Potential.contaminant != "+")
#' example_unique <- unique_names(example, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(example_unique))
#' example_exprset <- into_exprset(example_unique, columns)
#'
#' example_vsn <- norm_vsn(example_exprset)
#' example_filter <- miss_val_filter(example_vsn, thr = 0)
#' example_imp <- imputation_perseus(example_filter, shift = 1.8, scale = 0.3)
#' @export
imputation_perseus <- function(data, shift = 1.8, scale = 0.3) {
  stat <- exprs(data) %>% data.frame() %>% mutate(name = rownames(.)) %>% gather(samples, value, 1:(ncol(.)-1)) %>% filter(!is.na(value))  %>% group_by(samples) %>%
    summarise(mean = mean(value), median = median(value), sd = sd(value), n = n(), infin = nrow(exprs(data))-n)
  for(a in 1:nrow(stat)) {
    exprs(data)[is.na(exprs(data)[,stat$samples[a]]),stat$samples[a]] <- rnorm(stat$infin[a] , mean = stat$median[a] - shift * stat$sd[a], sd = stat$sd[a] * scale)
  }
  fData(data)$mean <- rowMeans(exprs(data))
  return(data)
}

#' Impute missing values
#'
#' \code{imputation_MSn} imputes missing values based on \code{\link{impute}}.
#'
#' @param data ExpressionSet, Data object for which missing values will be imputed.
#' @param fun Character, Function used for data imputation based on \code{\link{impute}}.
#' @return An imputed ExpressionSet object.
#' @examples
#' example <- UbIA_MS %>% filter(Reverse != "+", Potential.contaminant != "+")
#' example_unique <- unique_names(example, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(example_unique))
#' example_exprset <- into_exprset(example_unique, columns)
#'
#' example_vsn <- norm_vsn(example_exprset)
#' example_filter <- miss_val_filter(example_vsn, thr = 0)
#' example_impute_leftShifted <- imputation_MSn(example_filter, fun = "QRILC")
#' example_impute_nearestNeighbors <- imputation_MSn(example_filter, fun = "knn")
#' @export
imputation_MSn <- function(data, fun) {
  MSnSet_data <- as(data, "MSnSet")
  MSnSet_imputed <- impute(MSnSet_data, method = fun)
  exprs(data) <- exprs(MSnSet_imputed)
  fData(data)$mean <- rowMeans(exprs(data))
  return(data)
}

#' ANOVA and post-hoc Tukey
#'
#' \code{anova_tukey} performs an Analysis of Variance fit (\code{\link{aov}}) and subsequent Tukey Honest Significant Differences analysis (\code{\link{TukeyHSD}}).
#'
#' @param data ExpressionSet, Data object for which the variance will be analyzed.
#' @param control Character, The sample name to which the contrasts are generated (the control sample would be most appropriate).
#' @param type all" or "control" The type of contrasts that will be generated.
#' @return An ExpressionSet object containing FDR estimates of differential expression.
#' @examples
#' example <- UbIA_MS %>% filter(Reverse != "+", Potential.contaminant != "+")
#' example_unique <- unique_names(example, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(example_unique))
#' example_exprset <- into_exprset(example_unique, columns)
#'
#' example_vsn <- norm_vsn(example_exprset)
#' example_filter <- miss_val_filter(example_vsn, thr = 0)
#' example_impute <- imputation_MSn(example_filter, fun = "QRILC")
#'
#' example_anova_tukey <- anova_tukey(example_impute, "Con_", "control")
#' @export
anova_tukey <- function(data, control, type) {
  long <- exprs(data) %>% data.frame() %>% mutate(name = rownames(.)) %>% gather(ID, val, 1:(ncol(.)-1)) %>% left_join(., pData(data))

  cat("  ANOVA test \n")
  anova_p <- long %>% group_by(name) %>% do(anova = aov(val ~ sample, data = .) %>% summary(.) %>% .[[1]] %>% .$`Pr(>F)` %>% .[1])
  anova_p %<>% ungroup(name) %>% mutate(p = unlist(anova), padj = p.adjust(p, method = "BH")) %>% select(-anova) %>% data.frame()
  fData(data) <- left_join(fData(data), anova_p)

  cat("\n  Post-hoc test (Tukey) \n")
  tukey <- long %>% group_by(name) %>% do(tukey = aov(val ~ sample, data = .) %>% TukeyHSD(.) %>% .$sample %>% .[,c(1,4)])
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
  fData(data) <- merge(fData(data), final, by="name")
  rownames(fData(data)) <- fData(data)$name
  return(data)
}

#' Linear model test
#'
#' \code{linear_model} performs a differential expression test based on linear models and empherical Bayes statistics (\code{\link{limma}}).
#'
#' @param data ExpressionSet, Data object for which the variance will be analyzed.
#' @param control Character, The sample name to which the contrasts are generated (the control sample would be most appropriate).
#' @param type all" or "control" The type of contrasts that will be generated.
#' @return An ExpressionSet object containing FDR estimates of differential expression.
#' @examples
#' example <- UbIA_MS %>% filter(Reverse != "+", Potential.contaminant != "+")
#' example_unique <- unique_names(example, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(example_unique))
#' example_exprset <- into_exprset(example_unique, columns)
#'
#' example_vsn <- norm_vsn(example_exprset)
#' example_filter <- miss_val_filter(example_vsn, thr = 0)
#' example_impute <- imputation_MSn(example_filter, fun = "QRILC")
#'
#' example_lm <- linear_model(example_impute, "Con_", "control")
#' @export
linear_model <- function(data, control, type) {
  samples <- factor(pData(data)$sample)
  design <- model.matrix(~ 0 + samples)
  colnames(design) <- gsub("samples", "", colnames(design))

  if(type == "all") {
    cntrst <- apply(combn(unique(pData(data)$sample), 2), 2, function(x) paste(x, collapse = " - "))
    flip <- grep(paste("^", control, sep =""), cntrst)
    if(length(flip) >= 1) {
      cntrst[flip] %<>% gsub(paste(control, "- ", sep = " "), "", .) %>% paste(" - ", control, sep = "")
    }
  }
  if(type == "control") {
    cntrst <- paste(colnames(design)[!colnames(design) %in% control], control, sep = " - ")
  }
  eB_fit <- eBayes(contrasts.fit(lmFit(data, design = design), makeContrasts(contrasts = cntrst, levels = design)))

  retrieve_fun <- function(comp, fit = eB_fit){
    res <- topTable(fit, sort.by = "t", coef = comp, number = Inf)
    fdr_res <- fdrtool(res$t, plot = FALSE, verbose = FALSE)
    res$qval <- fdr_res$qval
    res$lfdr <- fdr_res$lfdr
    res$comparison <- rep(comp, dim(res)[1])
    return(res)
  }

  limma_res <- map_df(cntrst, retrieve_fun)

  limma_res_small <- limma_res %>% select(name, logFC, qval, comparison) %>% mutate(comparison = gsub(" - ", "-", comparison))
  table_diff <- limma_res_small %>% select(name, logFC, comparison) %>% spread(comparison, logFC)
  colnames(table_diff)[2:ncol(table_diff)] <- paste(colnames(table_diff)[2:ncol(table_diff)], "diff", sep = "_")
  table_padj <- limma_res_small %>% select(name, qval, comparison) %>% spread(comparison, qval)
  colnames(table_padj)[2:ncol(table_padj)] <- paste(colnames(table_padj)[2:ncol(table_padj)], "p.adj", sep = "_")
  table <- left_join(table_diff, table_padj, by = "name")
  fData(data) <- merge(fData(data), table, by = "name")
  rownames(fData(data)) <- fData(data)$name
  return(data)
}


#' Denote significant proteins
#'
#' \code{cutoffs} denotes significant proteins based on defined cutoffs.
#'
#' @param data ExpressionSet, Data object which will be filtered for significant proteins.
#' @param alpha Integer, sets the threshold for the false discovery rate (FDR).
#' @param lfc Integer, sets the threshold for the log fold change (lfc).
#' @return An ExpressionSet object containing a variable "sign" denoting the significant proteins.
#' @examples
#' example <- UbIA_MS %>% filter(Reverse != "+", Potential.contaminant != "+")
#' example_unique <- unique_names(example, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(example_unique))
#' example_exprset <- into_exprset(example_unique, columns)
#' example_vsn <- norm_vsn(example_exprset)
#' example_filter <- miss_val_filter(example_vsn, thr = 0)
#' example_impute <- imputation_MSn(example_filter, fun = "QRILC")
#'
#' example_lm <- linear_model(example_impute, "Con_", "control")
#' example_sign <- cutoffs(example_lm, alpha = 0.05, lfc = 1)
#' significant_proteins <- example_sign[fData(example_sign)$sign == "+"]
#' nrow(significant_proteins)
#' @export
cutoffs <- function(data, alpha = 0.05, lfc = 1) {
  feat_data <- fData(data)
  cols_p <- grep("_p.adj",colnames(feat_data))
  cols_diff <- grep("_diff", colnames(feat_data))
  if(length(cols_p) > 1) {
    sign_df <- ifelse(feat_data[,cols_p] %>% apply(., 2, function(x) x <= alpha) & feat_data[,cols_diff] %>% apply(., 2, function(x) x >= lfc | x <= -lfc), "+", "") %>% data.frame()
    sign_df %<>% mutate(sign = ifelse(apply(., 1, function(x) any(x == "+")),"+",""))
  } else {
    sign_df <- ifelse(feat_data[,cols_p] <= alpha & feat_data[,cols_diff] >= lfc | feat_data[,cols_p] <= alpha & feat_data[,cols_diff] <= -lfc, "+", "") %>% data.frame()
    colnames(sign_df) <- "sign"
  }
  colnames(sign_df) %<>% gsub("_p.adj","_sign",.)
  fData(data) <- cbind(fData(data), sign_df)
  return(data)
}

#' Generate a results table
#'
#' \code{results} generates a results table (data.frame) from a ExpressionSet object which has been generated by \code{\link{linear_model}} followed by \code{\link{cutoffs}}.
#'
#' @param data ExpressionSet, Data object which has been generated by \code{\link{linear_model}} and \code{\link{cutoffs}}.
#' @return An data.frame object containing all results variables from the performed analysis.
#' @examples
#' example <- UbIA_MS %>% filter(Reverse != "+", Potential.contaminant != "+")
#' example_unique <- unique_names(example, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(example_unique))
#' example_exprset <- into_exprset(example_unique, columns)
#' example_vsn <- norm_vsn(example_exprset)
#' example_filter <- miss_val_filter(example_vsn, thr = 0)
#' example_impute <- imputation_MSn(example_filter, fun = "QRILC")
#'
#' example_lm <- linear_model(example_impute, "Con_", "control")
#' example_sign <- cutoffs(example_lm, alpha = 0.05, lfc = 1)
#' example_results <- results(example_sign)
#' glimpse(example_results)
#' @export
results <- function(data) {
  feat_data <- fData(data)
  centered <- exprs(data) - fData(data)$mean
  centered %<>%  data.frame(.) %>% rownames_to_column(.) %>% gather(ID, val, 2:ncol(.)) %>% left_join(., pData(data), by = "ID")
  centered %<>% group_by(rowname, sample) %>% summarize(val = mean(val)) %>% mutate(val = signif(val, digits = 3)) %>% spread(sample, val)
  colnames(centered)[2:ncol(centered)] %<>%  paste(., "_centered", sep = "")

  ratio <- feat_data[,grep("_diff", colnames(feat_data))] %>% signif(., digits = 3) %>% rownames_to_column(.)
  colnames(ratio)[2:ncol(ratio)] %<>% gsub("_diff", "_ratio", .)
  df <- left_join(ratio, centered, by = "rowname")

  ids <- feat_data %>% select(name, ID)

  pval <- feat_data[,grep("p.adj|sign", colnames(feat_data))] %>% rownames_to_column(.)
  pval[,grep("p.adj", colnames(pval))] %<>%  signif(., digits = 3) %>% format(., scientific = T)

  table <- left_join(ids, pval, by = c("name" = "rowname"))
  table <- left_join(table, df, by = c("name" = "rowname"))
  return(table)
}
