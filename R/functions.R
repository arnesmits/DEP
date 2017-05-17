#' Make unique names
#'
#' \code{make_unique} generates unique identifiers
#' for a proteomics dataset based on "name" and "id" columns.
#'
#' @param data Data.frame,
#' Protein table for which unique names will be created.
#' @param name Character,
#' Name of the column containing feature names.
#' @param ids Character,
#' Name of the column containing feature IDs.
#' @param delim Character,
#' Delimiter separating the feature names within on protein group.
#' @return A data.frame with the additional variables
#' "name" and "ID" containing unique names and identifiers, respectively.
#' @examples
#' data <- UbiLength
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#' @export
make_unique <- function(data, name, ids, delim = ";") {
  # Show error if inputs are not the required classes
  assertthat::assert_that(is.data.frame(data),
                          is.character(name),
                          is.character(ids),
                          is.character(delim))

  # Show error if inputs do not contain required columns
  if(length(grep(paste("^", name, "$", sep = ""), colnames(data))) < 1) {
    stop(paste0("'", name, "' is not a column in '",
                deparse(substitute(data)), "'."),
         call. = FALSE)
  }
  if(length(grep(paste("^", ids, "$", sep = ""), colnames(data))) < 1) {
    stop(paste0("'", ids, "' is not a column in '",
                deparse(substitute(data)), "'."),
         call. = FALSE)
  }

  # If input is a tibble, convert to data.frame
  if(tibble::is.tibble(data)) data <- as.data.frame(data)

  # Select the name and id columns,
  # take the first identifier per row and make unique names.
  # If there is no name, the ID will be taken.
  names <- data %>% select(matches(paste("^", name, "$", sep = "")),
                           matches(paste("^", ids, "$", sep = ""))) %>%
    mutate(name = gsub(paste(delim, ".*", sep = ""), "", .[, 1]),
           ID = gsub(paste(delim, ".*", sep = ""), "", .[, 2]),
           name = make.unique(ifelse(name == "" | is.na(name), ID, name)))
  data_unique <- left_join(data, names)
  return(data_unique)
}

#' Data.frame to SummarizedExperiment object
#' conversion using an experimental design
#'
#' \code{make_se} creates a SummarizedExperiment object
#' based on two data.frames: the protein table and experimental design.
#'
#' @param data Data.frame,
#' Protein table with unique names annotated in the 'name' column.
#' @param columns Integer vector,
#' Column numbers indicating the columns containing the assay data.
#' @param expdesign Data.frame,
#' Experimental design with 'label', 'condition'
#' and 'replicate' information.
#' See \code{\link{UbiLength_ExpDesign}} for an example experimental design.
#' @return A SummarizedExperiment object
#' with log2-transformed values.
#' @examples
#' data <- UbiLength
#' data <- data[data$Reverse != "+" & data$Potential.contaminant != "+",]
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(data_unique))
#' exp_design <- UbiLength_ExpDesign
#' se <- make_se(data_unique, columns, exp_design)
#' @export
make_se <- function(data, columns, expdesign) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(is.data.frame(data),
                          is.integer(columns),
                          is.data.frame(expdesign))

  # Show error if inputs do not contain required columns
  if(any(!c("name", "ID") %in% colnames(data))) {
    stop(paste0("'name' and/or 'ID' columns are not present in '",
                deparse(substitute(data)),
                "'.\nRun make_unique() to obtain the required columns."),
         call. = FALSE)
  }
  if(any(!c("label", "condition", "replicate") %in% colnames(expdesign))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in the experimental design.",
         call. = FALSE)
  }
  if(any(!apply(data[, columns], 2, is.numeric))) {
    stop(paste0("specified columns (", columns,
                ") should be numeric.\nRun make_se_parse() with the appropriate columns as argument."),
         call. = FALSE)
  }

  # If input is a tibble, convert to data.frame
  if(tibble::is.tibble(data)) data <- as.data.frame(data)
  if(tibble::is.tibble(expdesign)) expdesign <- as.data.frame(expdesign)

  # Select the assay data
  rownames(data) <- data$name
  raw <- data[, columns]
  raw[raw == 0] <- NA
  raw <- log2(raw)

  # Generate the colData from the experimental design
  # and match these with the assay data
  expdesign <- mutate(expdesign, condition = make.names(condition)) %>%
    unite(ID, condition, replicate, remove = FALSE)
  rownames(expdesign) <- expdesign$ID
  colnames(raw)[lapply(expdesign$label,
                       function(x) grep(x, colnames(raw)))
                %>% unlist()] <- expdesign$ID
  raw <- raw[, !is.na(colnames(raw))][rownames(expdesign)]

  # Select the rowData
  row_data <- data[, -columns]
  rownames(row_data) <- row_data$name

  # Generate the SummarizedExperiment object
  dataset <- SummarizedExperiment(assays = as.matrix(raw),
                                  colData = expdesign,
                                  rowData = row_data)
  return(dataset)
}

#' Obtain the longest common prefix
#'
#' \code{get_prefix} returns the longest common prefix
#' of the supplied words.
#'
#' @param words Character vector,
#' A list of words.
#' @return A character vector containing the prefix.
#' @examples
#' data <- UbiLength
#' columns <- grep("LFQ.", colnames(data))
#'
#' names <- colnames(data[, columns])
#' get_prefix(names)
#' @export
get_prefix <- function(words) {
  # Show error if input is not the required class
  assertthat::assert_that(is.character(words))

  # Show error if 'words' contains 1 or less elements
  if(length(words) <= 1) {
    stop(paste0("'", deparse(substitute(words)),
                "' should contain more than one element."))
  }

  # Show error if 'words' contains NA
  if(any(is.na(words))) {
    stop(paste0("'", deparse(substitute(words)),
                "' contains NA(s)."))
  }

  # Truncate words to smallest name
  minlen <- min(nchar(words))
  truncated <- substr(words, 1, minlen)

  # Show error if one of the elements is shorter than one character
  if(minlen <= 1) {
    stop("At least one of the elements is too short.")
  }

  # Get identifical characters
  mat <- data.frame(strsplit(truncated, ""), stringsAsFactors = FALSE)
  identical <- apply(mat, 1, function(x) length(unique(x)) == 1)

  # Obtain the longest common prefix
  prefix <- as.logical(cumprod(identical))
  paste(mat[prefix, 1], collapse = "")
}

#' Data.frame to SummarizedExperiment object
#' conversion using parsing from column names
#'
#' \code{make_se_parse} creates a SummarizedExperiment object
#' based on a single data.frame.
#'
#' @param data Data.frame,
#' Protein table with unique names annotated in the 'name' column.
#' @param columns Integer vector,
#' Column numbers indicating the columns containing the assay data.
#' @return A SummarizedExperiment object
#' with log2-transformed values.
#' @examples
#' data <- UbiLength
#' data <- data[data$Reverse != "+" & data$Potential.contaminant != "+",]
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(data_unique))
#' se <- make_se_parse(data_unique, columns)
#' @export
make_se_parse <- function(data, columns) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(is.data.frame(data), is.integer(columns))

  # Show error if inputs do not contain required columns
  if(any(!c("name", "ID") %in% colnames(data))) {
    stop(paste0("'name' and/or 'ID' columns are not present in '",
                deparse(substitute(data)),
                "'.\nRun make_unique() to obtain the required columns."),
         call. = FALSE)
  }
  if(any(!apply(data[, columns], 2, is.numeric))) {
    stop(paste0("specified columns (", columns, ") should be numeric.
                \nRun make_se_parse() with the appropriate columns as argument."),
         call. = FALSE)
  }

  # If input is a tibble, convert to data.frame
  if(tibble::is.tibble(data)) data <- as.data.frame(data)

  # Select the assay values
  rownames(data) <- data$name
  raw <- data[, columns]
  raw[raw == 0] <- NA
  raw <- log2(raw)
  colnames(raw) <- gsub(get_prefix(colnames(raw)), "",
                        colnames(raw)) %>% make.names()

  # Select the rowData
  row_data <- data[, -columns]
  rownames(row_data) <- row_data$name

  # Generate the colData
  col_data <- data.frame(label = colnames(raw), stringsAsFactors = FALSE) %>%
    mutate(condition = substr(label, 1, nchar(label) - 1),
           replicate = substr(label, nchar(label), nchar(label))) %>%
    unite(ID, condition, replicate, remove = FALSE)
  rownames(col_data) <- col_data$ID
  colnames(raw) <- col_data$ID[
    lapply(col_data$label, function(x) grep(x, colnames(raw)))
    %>% unlist()]
  raw <- raw[, !is.na(colnames(raw))]

  # Generate the SummarizedExperiment object
  dataset <- SummarizedExperiment(assays = as.matrix(raw),
                                  colData = col_data,
                                  rowData = row_data)
  return(dataset)
}

#' Filter on missing values
#'
#' \code{filter_missval} filters a proteomics dataset based on missing values.
#'
#' @param data SummarizedExperiment,
#' Proteomics data with unique names and identifiers
#' annotated in 'name' and 'ID' columns.
#' The appropriate columns and objects can be generated
#' using \code{\link{make_se}} or \code{\link{make_se_parse}}.
#' @param thr Integer,
#' Sets the threshold for the allowed number of missing values per condition.
#' @return A filtered SummarizedExperiment object.
#' @examples
#' data <- UbiLength
#' data <- data[data$Reverse != "+" & data$Potential.contaminant != "+",]
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(data_unique))
#' exp_design <- UbiLength_ExpDesign
#' se <- make_se(data_unique, columns, exp_design)
#'
#' stringent_filter <- filter_missval(se, thr = 0)
#' less_stringent_filter <- filter_missval(se, thr = 1)
#' @export
filter_missval <- function(data, thr = 0) {
  # Show error if inputs are not the required classes
  if(is.integer(thr)) thr <- as.numeric(thr)
  assertthat::assert_that(inherits(data, "SummarizedExperiment"),
                          is.numeric(thr))

  # Show error if inputs do not contain required columns
  if(any(!c("name", "ID") %in% colnames(rowData(data)))) {
    stop(paste0("'name' and/or 'ID' columns are not present in '",
                deparse(substitute(data)),
                "'.\nRun make_unique() and make_se() to obtain the required columns."),
         call. = FALSE)
  }
  if(any(!c("label", "condition", "replicate") %in% colnames(colData(data)))) {
    stop(paste0("'label', 'condition' and/or 'replicate' columns are not present in '",
                deparse(substitute(data)),
                "'.\nRun make_se() or make_se_parse() to obtain the required columns."),
         call. = FALSE)
  }
  if(thr < 0 | thr > max(colData(data)$replicate)) {
    stop("invalid filter threshold applied;\nRun filter_missval() with a threshold ranging from 0 to the number of replicates")
  }

  # Make assay values binary (1 = valid value)
  bin_data <- assay(data)
  bin_data[!is.na(assay(data))] <- 1
  bin_data[is.na(assay(data))] <- 0

  # Filter data on the maximum allowed number of
  # missing values per condition (defined by thr)
  keep <- bin_data %>%
    data.frame() %>%
    rownames_to_column(.) %>%
    gather(ID, value, 2:ncol(.)) %>%
    left_join(., data.frame(colData(data)), by = "ID") %>%
    group_by(rowname, condition) %>%
    summarize(miss_val = n() - sum(value)) %>%
    filter(miss_val <= thr) %>%
    spread(condition, miss_val)
  filt <- data[keep$rowname, ]
  return(filt)
}

#' Normalization using vsn
#'
#' \code{normalize_vsn} performs variance stabilizing transformation
#' using the \code{\link[vsn]{vsn-package}}.
#'
#' @param data SummarizedExperiment,
#' Proteomics data with log2-transformed values
#' (as data obtained from \code{\link{make_se}}.
#' @return A normalized SummarizedExperiment object.
#' @examples
#' data <- UbiLength
#' data <- data[data$Reverse != "+" & data$Potential.contaminant != "+",]
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(data_unique))
#' exp_design <- UbiLength_ExpDesign
#' se <- make_se(data_unique, columns, exp_design)
#'
#' filt <- filter_missval(se, thr = 0)
#' norm <- normalize_vsn(filt)
#' @export
normalize_vsn <- function(data) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(data, "SummarizedExperiment"))

  # Variance stabilization transformation on assay data
  data_vsn <- data
  vsn.fit <- vsn::vsnMatrix(2 ^ assay(data_vsn))
  assay(data_vsn) <- vsn::predict(vsn.fit, 2 ^ assay(data_vsn))
  return(data_vsn)
}

#' Imputation by random draws from a manually defined distribution
#'
#' \code{manual_impute} imputes missing values in a proteomics dataset
#' by random draws from a manually defined distribution.
#'
#' @param data SummarizedExperiment,
#' Proteomics data for which missing values will be imputed.
#' @param shift Numeric,
#' Sets the left-shift of the distribution in standard deviation from
#' the mean of the original distribution.
#' @param scale Numeric,
#' Sets the width of the distribution relative to the
#' standard deviation of the original distribution.
#' @return An imputed SummarizedExperiment object.
#' @examples
#' data <- UbiLength
#' data <- data[data$Reverse != "+" & data$Potential.contaminant != "+",]
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(data_unique))
#' exp_design <- UbiLength_ExpDesign
#' se <- make_se(data_unique, columns, exp_design)
#'
#' filt <- filter_missval(se, thr = 0)
#' norm <- normalize_vsn(filt)
#'
#' imputed_manual <- impute(norm, fun = "man", shift = 1.8, scale = 0.3)
#' @export
manual_impute <- function(data, scale = 0.3, shift = 1.8) {
  if(is.integer(scale)) scale <- is.numeric(scale)
  if(is.integer(shift)) shift <- is.numeric(shift)
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(data, "SummarizedExperiment"),
                          is.numeric(scale),
                          is.numeric(shift))

  # Get descriptive parameters of the current sample distributions
  stat <- assay(data) %>%
    data.frame() %>%
    rownames_to_column(.) %>%
    gather(samples, value, 2:ncol(.)) %>%
    filter(!is.na(value))  %>% group_by(samples) %>%
    summarise(mean = mean(value),
              median = median(value),
              sd = sd(value),
              n = n(),
              infin = nrow(assay(data)) - n)
  # Impute missing values by random draws from a distribution
  # which is left-shifted by parameters 'shift' * sd and scaled by parameter 'scale' * sd.
  for (a in 1:nrow(stat)) {
    assay(data)[is.na(assay(data)[, stat$samples[a]]), stat$samples[a]] <-
      rnorm(stat$infin[a],
            mean = stat$median[a] - shift * stat$sd[a],
            sd = stat$sd[a] * scale)
  }
  return(data)
}

#' SummarizedExperiment to MSnSet object conversion
#'
#' \code{se2msn} generates a MSnSet object from a SummarizedExperiment object.
#'
#' @param data SummarizedExperiment,
#' Object which will be turned into a MSnSet object.
#' @return A MSnSet object.
#' @examples
#' data <- UbiLength
#' data <- data[data$Reverse != "+" & data$Potential.contaminant != "+",]
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(data_unique))
#' exp_design <- UbiLength_ExpDesign
#' se <- make_se(data_unique, columns, exp_design)
#'
#' data_msn <- se2msn(se)
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
#' \code{impute} imputes missing values in a proteomics dataset.
#'
#' @param data SummarizedExperiment,
#' Proteomics data with unique names and identifiers
#' annotated in 'name' and 'ID' columns.
#' The appropriate columns and objects can be generated
#' using \code{\link{make_se}} or \code{\link{make_se_parse}}.
#' @param fun "man", "bpca", "knn", "QRILC", "MLE", "MinDet",
#' "MinProb", "min", "zero", "mixed" or "nbavg",
#' Function used for data imputation based on \code{\link{manual_impute}}
#' and \code{\link[MSnbase]{impute}}.
#' @param ... Additional arguments for imputation functions as depicted in
#' \code{\link{manual_impute}} and \code{\link[MSnbase]{impute}}.
#' @return An imputed SummarizedExperiment object.
#' @examples
#' data <- UbiLength
#' data <- data[data$Reverse != "+" & data$Potential.contaminant != "+",]
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(data_unique))
#' exp_design <- UbiLength_ExpDesign
#' se <- make_se(data_unique, columns, exp_design)
#'
#' filt <- filter_missval(se, thr = 0)
#' norm <- normalize_vsn(filt)
#'
#' imputed_MinProb <- impute(norm, fun = "MinProb", q = 0.05)
#' imputed_QRILC <- impute(norm, fun = "QRILC")
#'
#' imputed_knn <- impute(norm, fun = "knn", k = 10, rowmax = 0.9)
#' imputed_MLE <- impute(norm, fun = "MLE")
#'
#' imputed_manual <- impute(norm, fun = "man", shift = 1.8, scale = 0.3)
#' @export
impute <- function(data, fun, ...) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(data, "SummarizedExperiment"),
                          is.character(fun))

  # Show error if inputs do not contain required columns
  if(any(!c("name", "ID") %in% colnames(rowData(data)))) {
    stop(paste0("'name' and/or 'ID' columns are not present in '",
                deparse(substitute(data)),
                "'.\nRun make_unique() and make_se() to obtain the required columns."), call. = FALSE)
  }
  if(!fun %in% c("man", MSnbase::imputeMethods())) {
    stop(paste("run imputation() with a valid function;\nValid functions are ",
               paste(c("man", MSnbase::imputeMethods()),
                     collapse = "', '"), "", sep = "'"))
  }

  # if the "man" function is selected, use the manual impution method
  if(fun == "man") {
    data <- manual_impute(data, ...)
  }
  # else use the MSnSet::impute function
  else {
    MSnSet_data <- se2msn(data)
    MSnSet_imputed <- MSnbase::impute(MSnSet_data, method = fun, ...)
    assay(data) <- MSnbase::exprs(MSnSet_imputed)
  }
  return(data)
}

#' Differential enrichment test
#'
#' \code{test_diff} performs a differential enrichment test based on
#' protein-wise linear models and empirical Bayes
#' statistics using \code{\link[limma]{limma}}.
#'
#' @param data SummarizedExperiment,
#' Proteomics data with unique names and identifiers
#' annotated in 'name' and 'ID' columns.
#' Additionally, the colData should contain sample annotation including
#' 'label', 'condition' and 'replicate' columns.
#' The appropriate columns and objects can be generated
#' using \code{\link{make_se}} or \code{\link{make_se_parse}}.
#' @param control Character,
#' The condition to which contrasts are generated
#' (a control condition would be most appropriate).
#' @param type "all", "control" or "manual",
#' The type of contrasts that will be tested.
#' This can be all possible pairwise comparisons ("all"),
#' limited to the comparisons versus the control ("control"), or
#' manually defined contrasts ("manual").
#' @param test Character,
#' The contrasts that will be tested if type = "manual".
#' These should be formatted as "SampleA_vs_SampleB" or
#' c("SampleA_vs_SampleC", "SampleB_vs_SampleC").
#' @return A SummarizedExperiment object
#' containing FDR estimates of differential expression.
#' @examples
#' data <- UbiLength
#' data <- data[data$Reverse != "+" & data$Potential.contaminant != "+",]
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(data_unique))
#' exp_design <- UbiLength_ExpDesign
#' se <- make_se(data_unique, columns, exp_design)
#'
#' filt <- filter_missval(se, thr = 0)
#' norm <- normalize_vsn(filt)
#' imputed <- impute(norm, fun = "MinProb", q = 0.01)
#'
#' diff <- test_diff(imputed, "Ctrl", "control")
#' diff <- test_diff(imputed, "Ctrl", "manual",
#'     test = c("Ubi4_vs_Ctrl", "Ubi6_vs_Ctrl"))
#' @export
test_diff <- function(data, control, type, test = NULL) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(data, "SummarizedExperiment"),
                          is.character(control),
                          is.character(type))

  # Show error if inputs do not contain required columns
  if(any(!c("name", "ID") %in% colnames(rowData(data)))) {
    stop(paste0("'name' and/or 'ID' columns are not present in '",
                deparse(substitute(data)),
                "'.\nRun make_unique() and make_se() to obtain the required columns."), call. = FALSE)
  }
  if(any(!c("label", "condition", "replicate") %in% colnames(colData(data)))) {
    stop(paste0("'label', 'condition' and/or 'replicate' columns are not present in '",
                deparse(substitute(data)),
                "'.\nRun make_se() or make_se_parse() to obtain the required columns."), call. = FALSE)
  }
  # Show error if inputs are not valid
  if(!type %in% c("all", "control", "manual")) {
    stop("run test_diff() with a valid type.\nValid types are: 'all', 'control'", call. = FALSE)
  }
  if(!control %in% unique(colData(data)$condition)) {
    stop("run test_diff() with a valid control",
         paste0(".\nValid controls are: '",
                paste0(unique(colData(data)$condition),
                       collapse = "', '"), "'"), call. = FALSE)
  }

  # Make an appropriate design matrix
  conditions <- factor(colData(data)$condition)
  design <- model.matrix(~ 0 + conditions)
  colnames(design) <- gsub("conditions", "", colnames(design))

  # Generate contrasts to be tested
  # Either make all possible combinations ("all") or
  # only the contrast versus the control sample ("control")
  if(type == "all") {
    cntrst <- apply(combn(colnames(design), 2), 2,
                    function(x) paste(x, collapse = " - "))
    # Make sure that contrast containing
    # the control sample have the control as denominator
    flip <- grep(paste("^", control, sep = ""), cntrst)
    if(length(flip) >= 1) {
      cntrst[flip] <- cntrst[flip] %>%
        gsub(paste(control, "- ", sep = " "), "", .) %>%
        paste(" - ", control, sep = "")
    }
  }
  if(type == "control") {
    cntrst <- paste(colnames(design)[!colnames(design) %in% control],
                    control,
                    sep = " - ")
  }
  if(type == "manual") {
    if(is.null(test)) {
      stop("run test_diff(type = 'manual') with a 'test' argument")
    }
    assertthat::assert_that(is.character(test))

    if(any(!unlist(strsplit(test, "_vs_")) %in% conditions)) {
      stop("run test_diff() with valid contrasts in 'test'",
           paste0(".\nValid contrasts should contain combinations of: '",
                  paste0(unique(colData(data)$condition),
                         collapse = "', '"),
                  "'.\nFor example '",
                  paste0(unique(colData(data)$condition)[1],
                         "_vs_",
                         unique(colData(data)$condition)[2])
                  , "'."), call. = FALSE)
    }

    cntrst <- gsub("_vs_", " - ", test)

  }
  # Print tested contrasts
  cat("Tested contrasts: \n")
  print(gsub(" - ", "_vs_", cntrst))

  # Test for differential expression by empirical Bayes moderation
  # of a linear model on the predefined contrasts
  eB_fit <- eBayes(
    contrasts.fit(
      lmFit(assay(data), design = design),
      makeContrasts(contrasts = cntrst, levels = design)
    )
  )

  # function to retrieve the results of
  # the differential expression test using 'fdrtool'
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
  limma_res_small <- limma_res %>%
    select(rowname, logFC, qval, comparison) %>%
    mutate(comparison = gsub(" - ", "_vs_", comparison))
  # Obtain a wide table with all log2 fold changes
  table_diff <- limma_res_small %>%
    select(rowname, logFC, comparison) %>%
    spread(comparison, logFC)
  colnames(table_diff)[2:ncol(table_diff)] <-
    paste(colnames(table_diff)[2:ncol(table_diff)], "diff", sep = "_")
  # Obtain a wide table with all q-values
  table_padj <- limma_res_small %>%
    select(rowname, qval, comparison) %>% spread(comparison, qval)
  colnames(table_padj)[2:ncol(table_padj)] <-
    paste(colnames(table_padj)[2:ncol(table_padj)], "p.adj", sep = "_")
  # Join the two tables with the rowData
  table <- left_join(table_diff, table_padj, by = "rowname")
  rowData(data) <- merge(rowData(data), table, by.x = "name", by.y = "rowname")
  return(data)
}

#' Mark significant proteins
#'
#' \code{add_rejections} marks significant proteins based on defined cutoffs.
#'
#' @param diff SummarizedExperiment,
#' Proteomics dataset on which differential enrichment analysis
#' has been performed by \code{\link{test_diff}}.
#' @param alpha Numeric,
#' Sets the threshold for the adjusted P value.
#' @param lfc Numeric,
#' Sets the threshold for the log2 fold change.
#' @return A SummarizedExperiment object
#' annotated with logical columns indicating significant proteins.
#' @examples
#' data <- UbiLength
#' data <- data[data$Reverse != "+" & data$Potential.contaminant != "+",]
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(data_unique))
#' exp_design <- UbiLength_ExpDesign
#' se <- make_se(data_unique, columns, exp_design)
#'
#' filt <- filter_missval(se, thr = 0)
#' norm <- normalize_vsn(filt)
#' imputed <- impute(norm, fun = "MinProb", q = 0.01)
#'
#' diff <- test_diff(imputed, "Ctrl", "control")
#' signif <- add_rejections(diff, alpha = 0.05, lfc = 1)
#' @export
add_rejections <- function(diff, alpha = 0.05, lfc = 1) {
  # Show error if inputs are not the required classes
  if(is.integer(alpha)) alpha <- as.numeric(alpha)
  if(is.integer(lfc)) lfc <- as.numeric(lfc)
  assertthat::assert_that(inherits(diff, "SummarizedExperiment"),
                          is.numeric(alpha),
                          is.numeric(lfc))

  # Show error if inputs do not contain required columns
  if(any(!c("name", "ID") %in% colnames(rowData(diff)))) {
    stop(paste0("'name' and/or 'ID' columns are not present in '",
                deparse(substitute(diff)),
                "';\nRun make_unique() and make_se() to obtain the required columns"),
         call. = FALSE)
  }
  if(length(grep("_p.adj|_diff", colnames(rowData(diff)))) < 1) {
    stop(paste0("'[contrast]_diff' and/or '[contrast]_p.adj' columns are not present in '",
                deparse(substitute(diff)),
                "';\nRun test_diff() to obtain the required columns"),
         call. = FALSE)
  }

  # get all columns with adjusted p-values and log2 fold changes
  row_data <- rowData(diff)
  cols_p <- grep("_p.adj", colnames(row_data))
  cols_diff <- grep("_diff", colnames(row_data))

  # Mark differential expressed proteins by
  # applying alpha and log2FC parameters per protein
  if(length(cols_p) == 1) {
    rowData(diff)$significant <-
      row_data[, cols_p] <= alpha & row_data[, cols_diff] >= lfc |
      row_data[, cols_p] <= alpha & row_data[, cols_diff] <= -lfc
    rowData(diff)$contrast_significant <- rowData(diff)$significant
    colnames(rowData(diff))[ncol(rowData(diff))] <-
      gsub("p.adj", "significant", colnames(row_data)[cols_p])
  }
  if(length(cols_p) > 1) {
    sign_df <- row_data[, cols_p] %>%
      apply(., 2, function(x) x <= alpha) &
      row_data[, cols_diff] %>% apply(., 2, function(x) x >= lfc | x <= -lfc)
    sign_df <- cbind(sign_df,
                     significant = apply(sign_df, 1, function(x) any(x)))
    colnames(sign_df) <- gsub("_p.adj", "_significant", colnames(sign_df))

    sign_df <- cbind(name = row_data$name, as.data.frame(sign_df))
    rowData(diff) <- merge(rowData(diff), sign_df, by = "name")
  }
  return(diff)
}

#' Generate a results table
#'
#' \code{get_results} generates a results table from a proteomics dataset
#' on which differential enrichment analysis was performed.
#'
#' @param data SummarizedExperiment,
#' Proteomics dataset on which differential enrichment proteins
#' are annotated by \code{\link{add_rejections}}.
#' @return A data.frame object
#' containing all results variables from the performed analysis.
#' @examples
#' data <- UbiLength
#' data <- data[data$Reverse != "+" & data$Potential.contaminant != "+",]
#' data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")
#'
#' columns <- grep("LFQ.", colnames(data_unique))
#' exp_design <- UbiLength_ExpDesign
#' se <- make_se(data_unique, columns, exp_design)
#'
#' filt <- filter_missval(se, thr = 0)
#' norm <- normalize_vsn(filt)
#' imputed <- impute(norm, fun = "MinProb", q = 0.01)
#'
#' diff <- test_diff(imputed, "Ctrl", "control")
#' signif <- add_rejections(diff, alpha = 0.05, lfc = 1)
#'
#' results <- get_results(signif)
#' colnames(results)
#'
#' significant_proteins <- results[results$significant,]
#' nrow(significant_proteins)
#' head(significant_proteins)
#' @export
get_results <- function(data) {
  # Show error if inputs are not the required classes
  assert_that(inherits(data, "SummarizedExperiment"))

  # Show error if inputs do not contain required columns
  if (any(!c("name", "ID") %in% colnames(rowData(data)))) {
    stop(paste0("'name' and/or 'ID' columns are not present in '",
                deparse(substitute(data)),
                "'.\nRun make_unique() and make_se() to obtain the required columns."),
         call. = FALSE)
  }
  if (length(grep("_p.adj|_diff", colnames(rowData(data)))) < 1) {
    stop(paste0("'[contrast]_diff' and/or '[contrast]_p.adj' columns are not present in '",
                deparse(substitute(data)),
                "'.\nRun test_diff() to obtain the required columns."),
         call. = FALSE)
  }

  row_data <- data.frame(rowData(data))

  # Obtain average protein-centered enrichment values per condition
  rowData(data)$mean <- rowMeans(assay(data))
  centered <- assay(data) - rowData(data)$mean
  centered <- data.frame(centered) %>%
    rownames_to_column(.) %>%
    gather(ID, val, 2:ncol(.)) %>%
    left_join(., data.frame(colData(data)), by = "ID")
  centered <- group_by(centered, rowname, condition) %>%
    summarize(val = mean(val)) %>%
    mutate(val = signif(val, digits = 3)) %>%
    spread(condition, val)
  colnames(centered)[2:ncol(centered)] <-
    paste(colnames(centered)[2:ncol(centered)], "_centered", sep = "")

  # Obtain average enrichments of conditions versus the control condition
  ratio <- row_data %>%
    column_to_rownames("name") %>%
    select(ends_with("diff")) %>%
    signif(., digits = 3) %>%
    rownames_to_column(.)
  colnames(ratio)[2:ncol(ratio)] <-
    gsub("_diff", "_ratio", colnames(ratio)[2:ncol(ratio)])
  df <- left_join(ratio, centered, by = "rowname")

  # Select the adjusted p-values and significance columns
  pval <- row_data %>%
    column_to_rownames("name") %>%
    select(ends_with("p.adj"), ends_with("significant")) %>%
    rownames_to_column(.)
  pval[, grep("p.adj", colnames(pval))] <-
    format(signif(pval[, grep("p.adj", colnames(pval))]
                  , digits = 3),
           scientific = TRUE)

  # Join into a results table
  ids <- row_data %>% select(name, ID)
  table <- left_join(ids, pval, by = c("name" = "rowname"))
  table <- left_join(table, df, by = c("name" = "rowname")) %>%
    arrange(desc(significant))
  return(table)
}
