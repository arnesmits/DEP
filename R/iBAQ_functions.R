#' Merge iBAQ intensities
#'
#' \code{merge_ibaq} generates a data.frame with merged iBAQ intensities
#' for proteins with shared peptides.
#'
#' @param protein_unique Data.frame,
#' Protein table with unique names annotated in the 'name' column.
#' @param peptides Data.frame,
#' Peptides table from MaxQuant ('peptides.txt').
#' @return A data.frame with iBAQ intensities per protein group.
#' Every protein group contains peptides unique to this group.
#' @examples
#' data <- GFPip
#' data_unique <- make_unique(data, 'Gene.names', 'Protein.IDs')
#' peptides <- GFPip_pep
#'
#' ibaq <- merge_ibaq(data_unique, peptides)
#' colnames(ibaq)
#' head(ibaq)
#' @export
merge_ibaq <- function(protein_unique, peptides) {
    # Show error if inputs are not the required classes
    assertthat::assert_that(is.data.frame(protein_unique),
                            is.data.frame(peptides))

    # Show error if inputs do not contain required columns
    if (any(!c("name", "ID") %in% colnames(protein_unique))) {
        stop(paste0("'name' and/or 'ID' columns are not present in '",
                    deparse(substitute(protein_unique)),
                    "'.\nRun make_unique() to obtain the required columns."),
             call. = FALSE)
    }
    if (length(grep("iBAQ.", colnames(protein_unique))) < 1) {
        stop(paste0("'iBAQ' columns are not present in '",
                    deparse(substitute(protein_unique)), "'."),
             call. = FALSE)
    }
    if (any(!c("Protein.group.IDs", "Unique..Groups.") %in% colnames(peptides))) {
        stop(paste0("'Protein.group.IDs' and/or 'Unique..Groups.' columns are not present in '",
                    deparse(substitute(peptides)), "'."),
             call. = FALSE)
    }

    # If input is a tibble, convert to data.frame
    if(tibble::is.tibble(protein_unique)) protein_unique <- as.data.frame(protein_unique)
    if(tibble::is.tibble(peptides)) peptides <- as.data.frame(peptides)

    # Filter for peptides not unique to a single protein group
    # and sort their protein group IDs
    shared_pep <- peptides %>% filter(Unique..Groups. == "no")
    # Function to split, sort and collapse IDs
    split_sort_collapse <- function(x) {
      vapply(x, function(y) {
        feat <- unlist(strsplit(y, ";"))
        sorted <- paste(sort(feat), collapse = ";")
      }, character(1))
    }

    # Split, sort and collapse protein group IDs and keep unique combinations
    sorted_pep <- shared_pep %>%
      select(Protein.group.IDs) %>%
      mutate(Protein.group.IDs = split_sort_collapse(Protein.group.IDs)) %>%
      unique()

    # Expand the protein groups from a single column to a matrix
    max <- vapply(sorted_pep$Protein.group.IDs,
                  function(x) length(unlist(strsplit(x, ";"))),
                  numeric(1)) %>% max()
    shared <- separate(sorted_pep,
                       Protein.group.IDs,
                       paste0("X", 1:max),
                       sep = ";", fill = "right")

    # Combine all protein group IDs that have overlapping peptides
    # set first IDs group
    x <- list(shared[1, !is.na(shared[1, ])])
    for (i in 2:nrow(shared)) {
      # get IDs of new row and check whether they match
      #with any of the previous IDs
      ids <- shared[i, !is.na(shared[i, ])]
      if (any(ids %in% unlist(x))) {
        # Get the IDs that match
        match <- ids[ids %in% unlist(x)]
        if (length(match) == 1) {
          # Grep the rows of the single match
          rows <- grep(paste("\\b", match, "\\b", sep = ""), x)
        }
        if (length(match) > 1) {
          # Grep the rows of all matches
          rows <- lapply(match,
                         function(y) {
                           grep(paste("\\b", y, "\\b", sep = ""), x) }) %>%
            unlist() %>% unique()
        }
        if (length(rows) == 1) {
          # Combine the single matched ID group with the current IDs
          x[[rows]] <- unique(unlist(c(unlist(x[[rows]]), ids)))
        }
        if (length(rows) > 1) {
          # Combine all matched ID groups and the current IDs in
          # the first matched ID group and remove the others
          x[[rows[1]]] <- unique(unlist(
            c(unlist(lapply(rows, function(y) x[[y]])), ids)))
          for (j in 2:length(rows)) {
            x[[rows[j] - j + 2]] <- NULL
          }
        }
      }
      if (!any(ids %in% unlist(x))) {
        # In case there is no match with previous IDs groups,
        # make a new IDs group
        x[[(length(x) + 1)]] <- ids
      }
    }

    # Function to convert the heterogenous list to a list of data.frames
    list2mat <- function(list) {
        if (is.character(list)) {
            dat <- t(as.data.frame(as.numeric(list)))
        } else {
            dat <- t(as.data.frame(as.numeric(list)))
        }
        colnames(dat) <- paste("X", 1:ncol(dat), sep = "")
        return(dat)
    }

    # Get all shared IDs
    shared_ids <- lapply(x, list2mat) %>% unlist(.)
    # Generate a single data.frame from the list of IDs groups
    rows2merge <- lapply(x, list2mat) %>%
      lapply(as.data.frame) %>%
      bind_rows()

    columns <- grep("iBAQ.", colnames(protein_unique))
    # Function to merge protein groups with shared peptides
    merge_sum <- function(rows) {
        sub <- protein_unique %>%
          filter(grepl(paste("^", rows, "$", collapse = "|", sep = ""), id))
        sub[1, columns] <- colSums(sub[, columns])
        sub$name[1] <- paste(sort(unique(sub$name)), collapse = ";")
        return(sub[1, ] %>% select(name, columns))
    }
    merged <- apply(rows2merge, 1, merge_sum) %>% bind_rows()

    # Function to count peptides for merged protein groups
    merge_pep <- function(rows) {
        sub <- peptides %>%
          filter(grepl(paste("^", rows, "$", collapse = "|", sep = ""),
                       Protein.group.IDs))
        return(nrow(sub))
    }
    peps <- apply(rows2merge, 1, merge_pep)
    merged$Peptides <- peps

    # Select all protein groups that only have peptides unique to this group
    data_unique <- protein_unique %>% filter(!id %in% shared_ids)

    # Generate the final list of all protein groups
    final <- rbind(data_unique %>%
                     select(name, columns, Peptides), merged)

    max <- lapply(final$name,
                  function(x) strsplit(x, ";")[[1]] %>% length()) %>%
      unlist() %>% max()
    names <- final %>%
      separate(., name, paste0("name_", 1:max),
               sep = ";", fill = "right", remove = FALSE)
    return(names)
}

#' Relative stoichiometry
#'
#' \code{get_stoichiometry} calculates the relative stoichiometries of
#' all differentially enriched proteins.
#'
#' @param dep SummarizedExperiment,
#' Proteomics dataset in which differential enriched proteins
#' are annotated by \code{\link{add_rejections}}.
#' @param ibaq Data.frame,
#' iBAQ table generated by \code{\link{merge_ibaq}}.
#' @param contrast Character,
#' The specific contrast for which to calculate the stoichiometries.
#' @param bait Character,
#' The name of the protein to which all other proteins will be scaled.
#' @param level Numerical,
#' The level to which the bait will be scaled.
#' @return A data.frame with relative stoichiometry values.
#' @examples
#' # load data
#' data <- GFPip
#' exp_design <- GFPip_ExpDesign
#' peptides <- GFPip_pep
#'
#' # Make unique names and filter for reverse and contaminant hits
#' data_unique <- make_unique(data, 'Gene.names', 'Protein.IDs')
#' data_filt <- data_unique[data_unique$Reverse != '+' & data_unique$Contaminant != '+',]
#'
#' # Make a SummarizedExperiment and perform filtering on missing values,
#' # variance stabilization and imputation on this object
#' cols <- grep('LFQ', colnames(data_filt))
#' se <- make_se(data_filt, cols, exp_design)
#' filt <- filter_missval(se)
#' norm <- normalize_vsn(filt)
#' imputed <- impute(norm, 'MinProb')
#'
#' # Test for differential enriched proteins
#' diff <- test_diff(imputed, 'WT', 'control')
#' final <- add_rejections(diff, alpha = 0.05, lfc = 4.5)
#'
#' # Merge iBAQ intensities of proteins that have shared peptides
#' ibaq <- merge_ibaq(data_unique, peptides)
#'
#' # Calculate relative stoichiometry versus 'Suz12' in the 'GFP_vs_WT' contrast
#' stoi <- get_stoichiometry(final, ibaq, contrast = 'GFP_vs_WT', bait = 'Suz12')
#'
#' @export
get_stoichiometry <- function(dep, ibaq, contrast, bait, level = 1) {
    # Show error if inputs are not the required classes
    if (is.integer(level))
        level <- as.numeric(level)
    assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                            is.data.frame(ibaq),
                            is.character(contrast),
                            is.character(bait),
                            is.numeric(level))

    # Get rowData
    row_data <- rowData(dep)

    # Show error if inputs do not contain required columns
    if (any(!c("name", "ID") %in% colnames(row_data))) {
        stop(paste0("'name' and/or 'ID' columns are not present in '",
                    deparse(substitute(dep)), "'."),
             call. = FALSE)
    }
    if (length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
        stop(paste0("'[contrast]_p.adj' and '[contrast]_diff' columns are not present in '",
                    deparse(substitute(dep)),
                    "'.\nRun test_diff() to obtain the required columns."),
             call. = FALSE)
    }
    if (length(grep("_significant", colnames(row_data))) < 1) {
        stop(paste0("[contrast]_significant' columns are not present in '",
                    deparse(substitute(dep)),
                    "'.\nRun get_rejections() to obtain the required columns."),
             call. = FALSE)
    }
    if (length(grep("iBAQ.", colnames(ibaq))) < 1) {
        stop(paste0("'iBAQ' columns are not present in '",
                    deparse(substitute(ibaq)), "'."),
             call. = FALSE)
    }
    if (length(grep("name_", colnames(ibaq))) < 1) {
        stop(paste0("merge information is not present in '",
                    deparse(substitute(ibaq)),
                    "'.\nRun merge_ibaq() to obtain the required 'name_' columns."),
             call. = FALSE)
    }

    # Show error if an unvalid contrast is given
    if (length(grep(paste(contrast, "_diff", sep = ""),
                    colnames(row_data))) == 0) {
        valid_cntrsts <- row_data %>%
          data.frame() %>%
          select(ends_with("_diff")) %>%
          colnames(.) %>%
          gsub("_diff", "", .)
        valid_cntrsts_msg <-
          paste0("Valid contrasts are: '",
                 paste0(valid_cntrsts, collapse = "', '"), "'")
        stop("Not a valid contrast, please run `plot_volcano()` with a valid contrast as argument\n",
             valid_cntrsts_msg,
             call. = FALSE)
    }

    # Get significant proteins for the defined contrast
    col_signif <- grep(paste(contrast, "_significant", sep = ""),
                       colnames(row_data))
    row_data <- row_data[row_data[, col_signif], ]

    # Select these significant proteins from the iBAQ table
    names <- row_data$name
    rows <- ibaq %>% select(starts_with("name_")) %>%
      apply(., 2,
            function(x) grep(paste("^", names, "$", sep = "", collapse = "|"), x)) %>%
      unlist() %>% unique()
    sub <- ibaq[rows, ] %>% select(name, starts_with("iBAQ."))

    # Generate an annotated and tidy table
    ip <- gsub("_vs_.*", "", contrast)
    control <- gsub(".*_vs_", "", contrast)
    ibaq_anno <- colData(dep) %>%
      data.frame() %>%
      mutate(sample = paste("iBAQ.", label, sep = ""))
    long <- sub %>%
      gather(sample, iBAQ, 2:(ncol(.) - 1)) %>%
      left_join(., ibaq_anno, by = "sample")

    # Calculate mean background signal per protein
    ctrl <- long %>%
      group_by(name, condition) %>%
      summarise(mean_ctrl = mean(iBAQ)) %>%
      filter(condition == control) %>%
      select(name, mean_ctrl)

    # Subtract background and generate a wide table
    stoi <- long %>%
      left_join(., ctrl, by = "name") %>%
      filter(condition == ip) %>%
      mutate(iBAQ = iBAQ - mean_ctrl)
    wide <- stoi %>%
      select(name, iBAQ, ID) %>%
      spread(ID, iBAQ)

    # Show error if an unvalid bait is given
    if (length(grep(paste("^", bait, "$", sep = ""), wide$name)) == 0) {
        valid_baits_msg <- paste0("Valid baits are: '",
                                  paste0(wide$name, collapse = "', '"),
                                  "'")
        stop("Not a valid bait, please run `stoichiometry()` with a valid bait as argument\n",
             valid_baits_msg,
             call. = FALSE)
    }

    # Get the bait iBAQ values
    num <- wide %>%
      filter(name == bait) %>%
      .[, 2:ncol(.)]
    num <- num / level

    # Calculate the relative stoichiometries versus bait protein
    df <- wide %>%
      mutate_at(2:ncol(.), funs(./num$.)) %>%
      gather(ID, iBAQ, 2:ncol(.)) %>%
      left_join(., ibaq_anno, by = "ID")

    # Generate a table with the mean ± sd stoichiometries
    final <- df %>%
      group_by(name, condition) %>%
      summarize(stoichiometry = mean(iBAQ), sd = sd(iBAQ)) %>%
      arrange(desc(stoichiometry)) %>%
      ungroup(name)
    return(final)
}

#' Plot relative stoichiometries
#'
#' \code{plot_stoichiometry} plots a barplot of the
#' relative stoichiometries.
#'
#' @param protein_stoichiometry Data.frame,
#' Stoichiometry table generated by \code{\link{get_stoichiometry}}.
#' @param thr Numerical,
#' The stoichiometry threshold above which proteins will be plotted.
#' @param max_y Numerical,
#' The maximum scale of the y-axis.
#' @return A barplot (generated by \code{\link[ggplot2]{ggplot}}).
#' @examples
#' # load data
#' data <- GFPip
#' exp_design <- GFPip_ExpDesign
#' peptides <- GFPip_pep
#'
#' # Make unique names and filter for reverse and contaminant hits
#' data_unique <- make_unique(data, 'Gene.names', 'Protein.IDs')
#' data_filt <- data_unique[data_unique$Reverse != '+' & data_unique$Contaminant != '+',]
#'
#' # Make a SummarizedExperiment and perform filtering on missing values,
#' # variance stabilization and imputation on this object
#' cols <- grep('LFQ', colnames(data_filt))
#' se <- make_se(data_filt, cols, exp_design)
#' filt <- filter_missval(se)
#' norm <- normalize_vsn(filt)
#' imputed <- impute(norm, 'MinProb')
#'
#' # Test for differential enriched proteins
#' diff <- test_diff(imputed, 'WT', 'control')
#' final <- add_rejections(diff, alpha = 0.05, lfc = 4.5)
#'
#' # Merge iBAQ intensities of proteins that have shared peptides
#' ibaq <- merge_ibaq(data_unique, peptides)
#'
#' # Calculate relative stoichiometry versus 'Suz12' in the 'GFP_vs_WT' contrast
#' stoi <- get_stoichiometry(final, ibaq, contrast = 'GFP_vs_WT', bait = 'Suz12')
#' plot_stoichiometry(stoi)
#'
#' @export
plot_stoichiometry <- function(protein_stoichiometry, thr = 0.01, max_y = NULL) {
    # Show error if inputs are not the required classes
    if (is.integer(thr)) thr <- as.numeric(thr)
    assertthat::assert_that(is.data.frame(protein_stoichiometry),
                            is.numeric(thr))
    if (!is.null(max_y)) {
        if (is.integer(max_y)) max_y <- as.numeric(max_y)
        assertthat::assert_that(is.numeric(max_y))
    }

    # Show error if data does not contain the required columns
    if (any(!c("name", "condition", "stoichiometry", "sd") %in% colnames(protein_stoichiometry))) {
        stop(paste0("'name', 'condition', 'stoichiometry' and/or 'sd' columns are not present in '",
                    deparse(substitute(protein_stoichiometry)),
                    "'.\nRun get_stoichiometry() to obtain required columns"),
            call. = FALSE)
    }

    # Obtain a table with stoichiometries
    df <- protein_stoichiometry %>%
      filter(stoichiometry >= thr) %>%
      mutate(name = ifelse(nchar(name) > 20,
                           paste(substr(name, 1, 20), "...", sep = ""),
                           name),
             ymin = stoichiometry - sd,
             ymax = stoichiometry + sd)
    df$name <- parse_factor(df$name, levels = unique(df$name))

    # Get the bait name
    bait <- df %>%
      filter(stoichiometry == 1) %>%
      .$name %>%
      unique()

    # Set the y-axis limit
    if (is.null(max_y)) {
      max <- max(df$ymax)
    } else {
      max <- max_y
    }

    # Plot the barplot with errorbars
    ggplot(df, aes(x = name, y = stoichiometry)) +
      geom_col(colour = "black", fill = "grey") +
      geom_errorbar(aes(ymax = ymax, ymin = ymin), width = 0.3) +
      labs(title = unique(df$condition), x = "",
           y = paste0("Stoichiometry (vs ", bait, ")")) +
      ylim(0, max) +
      theme_DEP2()
}

#' Plot iBAQ values versus log fold change
#'
#' \code{plot_ibaq} plots a scatter plot of the
#' iBAQ intensities versus the LFQ fold changes.
#'
#' @param dep SummarizedExperiment object,
#' Proteomics dataset on which differential enriched proteins
#' are annotated by \code{\link{add_rejections}}.
#' @param contrast Character,
#' The specific contrast to plot.
#' @param labelsize Integer,
#' Sets the size of name labels.
#' @return A scatter plot
#' @examples
#' # load data
#' data <- GFPip
#' exp_design <- GFPip_ExpDesign
#' peptides <- GFPip_pep
#'
#' # Make unique names and filter for reverse and contaminant hits
#' data_unique <- make_unique(data, 'Gene.names', 'Protein.IDs')
#' data_filt <- data_unique[data_unique$Reverse != '+' & data_unique$Contaminant != '+',]
#'
#' # Make a SummarizedExperiment and perform filtering on missing values,
#' # variance stabilization and imputation on this object
#' cols <- grep('LFQ', colnames(data_filt))
#' se <- make_se(data_filt, cols, exp_design)
#' filt <- filter_missval(se)
#' norm <- normalize_vsn(filt)
#' imputed <- impute(norm, 'MinProb')
#'
#' # Test for differential enriched proteins
#' diff <- test_diff(imputed, 'WT', 'control')
#' final <- add_rejections(diff, alpha = 0.05, lfc = 4.5)
#'
#' # Plot iBAQ vs LFQ plot
#' plot_ibaq(final, 'GFP_vs_WT', labelsize = 3)
#'
#' @export
plot_ibaq <- function(dep, contrast, labelsize = 3) {
    # Show error if inputs are not the required classes
    if (is.integer(labelsize)) labelsize <- as.numeric(labelsize)
    assertthat::assert_that(inherits(dep, "SummarizedExperiment"),
                            is.character(contrast),
                            is.numeric(labelsize))

    # Show error if inputs do not contain required columns
    if (any(!c("name", "ID") %in% colnames(rowData(dep)))) {
        stop(paste0("'name' and/or 'ID' columns are not present in '",
                    deparse(substitute(dep)), "'."),
             call. = FALSE)
    }
    if (length(grep("_p.adj|_diff", colnames(rowData(dep)))) < 1) {
        stop(paste0("'[contrast]_diff' and '[contrast]_p.adj' columns are not present in '",
                    deparse(substitute(dep)),
                    "'.\nRun test_diff() to obtain the required columns."),
             call. = FALSE)
    }
    if (length(grep("_significant", colnames(rowData(dep)))) < 1) {
        stop(paste0("'[contrast]_significant' columns are not present in '",
                    deparse(substitute(dep)),
                    "'.\nRun get_rejections() to obtain the required columns."),
             call. = FALSE)
    }
    if (length(grep("iBAQ.", colnames(rowData(dep)))) < 1) {
        stop(paste0("iBAQ columns are not present in '",
                    deparse(substitute(dep)), "'."),
             call. = FALSE)
    }

    # Get rowData
    row_data <- rowData(dep) %>% data.frame()

    if (length(grep(paste(contrast, "_diff", sep = ""),
                    colnames(row_data))) == 0) {
        valid_cntrsts <- row_data %>%
          data.frame() %>%
          select(ends_with("_diff")) %>%
          colnames(.) %>%
          gsub("_diff", "", .)
        valid_cntrsts_msg <- paste0("Valid contrasts are: '",
                 paste0(valid_cntrsts, collapse = "', '"), "'")
        stop("Not a valid contrast, please run `plot_ibaq()` with a valid contrast as argument\n",
             valid_cntrsts_msg,
             call. = FALSE)
    }

    # Generate an annotated and tidy table
    ip <- gsub("_vs_.*", "", contrast)
    control <- gsub(".*_vs_", "", contrast)
    ibaq_anno <- colData(dep) %>%
      data.frame() %>%
      mutate(sample = paste("iBAQ.", label, sep = ""))
    long <- row_data %>%
      gather(sample, iBAQ_value, starts_with("iBAQ.")) %>%
      select(name, sample, Peptides, iBAQ_value) %>%
      left_join(., ibaq_anno, by = "sample")

    # Subtract background level
    stoi <- long %>%
      group_by(name, condition) %>%
      summarise(mean = mean(iBAQ_value)) %>%
      spread(condition, mean)
    stoi$ibaq <- stoi[[ip]] - stoi[[control]]

    # Select the log2 fold change and adjusted p values columns
    # for the specific contrast
    col_diff <- grep(paste(contrast, "_diff", sep = ""),
                     colnames(row_data))
    col_signif <- grep(paste(contrast, "_significant", sep = ""),
                       colnames(row_data))
    row_data$lfc <- row_data[, col_diff]
    row_data$signif <- row_data[, col_signif]

    # Get the final table for plotting
    final <- row_data %>%
      select(name, lfc, Peptides, signif) %>%
      left_join(., stoi, by = "name")

    # Differentially enriched proteins
    extra <- final %>% filter(signif)

    # Generate the scatter plot
    ggplot(final, aes(x = lfc, y = log10(abs(ibaq)), size = Peptides)) +
      geom_vline(xintercept = 0) + geom_point(shape = 1, col = "grey") +
      geom_point(data = extra, shape = 1, col = "black") +
        scale_size_continuous(breaks = seq(0, 70, 10), range = c(1, 10)) +
      scale_color_manual(values = c("grey", "black")) +
      ggrepel::geom_text_repel(data = extra,
                               aes(label = name),
                               size = labelsize,
                               point.padding = unit(0.3, "lines")) +
      geom_text(data = data.frame(), aes(x = c(Inf, -Inf),
                                         y = c(-Inf, -Inf),
                                         hjust = c(1, 0),
                                         vjust = c(-1, -1),
                                         label = c(ip, control),
                                         fontface = "bold"),
                size = 5) +
      labs(title = contrast, x = "Fold change (log2)", y = "iBAQ intensity (log10)") +
      theme_DEP1()
}
